#!/usr/bin/env python
"""
Copyright 2018 

author:  <scott.wales@unimelb.edu.au>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""
from __future__ import print_function
import xarray
import subprocess
import numpy as np
import netCDF4
import pytest
import os
import glob
import sys

def run_nccopy(options, infiles, outdir):
    try:
        for f in infiles:
            fout = os.path.join(outdir,os.path.basename(f))
            s = subprocess.check_output(
                    ['nccopy'] + options + [f] + [fout],
                    stderr=subprocess.STDOUT)
            print(s.decode('utf-8'))
    except subprocess.CalledProcessError as e:
        print(e.stdout.decode('utf-8'))
        raise

# Run the collation program
def run_collate(inputs, output, np=2, args=[]):
    try:
        s = subprocess.check_output(
                ['mpirun', '-n', '%d'%np, './mppnccombine-fast', '--debug', '-o', str(output)] + inputs + args,
                stderr=subprocess.STDOUT)
        print(s.decode('utf-8'))
    except subprocess.CalledProcessError as e:
        print(e.stdout.decode('utf-8'))
        raise
    return xarray.open_dataset(str(output), engine='netcdf4', decode_times=False)

def split_dataset(data, split):
    out = []

    for start in range(0, len(data['x']), split['x']):
        d = data.isel(**{'x': slice(start, start+split['x'])})
        d['x'].attrs['domain_decomposition'] = [1,len(data['x']), 1+start, min(start+split['x'], len(data['x']))]

        out.append(d)

    return out


def split_file(tmpdir, data, split):
    i = 0
    infiles = []
    for start in range(0, len(data['x']), split['x']):
        infiles.append(str(tmpdir.join('in.%04d.nc'%i)))

        d = data.isel(**{'x': slice(start, start+split['x'])})
        d['x'].attrs['domain_decomposition'] = [1,len(data['x']), 1+start, min(start+split['x'], len(data['x']))]

        d.to_netcdf(infiles[-1],
                encoding = {'a': {
                    'chunksizes': (2,),
                    'zlib': True,
                    'shuffle': True,
                    'complevel': 4,
                    },
                    }
                )
        d.close()
        i+=1

    return infiles


def test_simple(tmpdir):
    inpath = str(tmpdir.join('test.nc'))
    outpath = tmpdir.join('out.nc')
    d = xarray.Dataset(
            {
                'a': (['x'], np.random.rand(5))
            },
            coords = {
                'x': np.arange(5),
            })

    d.to_netcdf(inpath, encoding={
            'a': {
                'chunksizes': (2,),
                'zlib': True,
                'shuffle': True,
                'complevel': 4,
                },
        })

    c = run_collate([inpath], outpath)

    assert (c.a.data == d.a.data).all()


def test_split_on_boundary(tmpdir):
    d = xarray.Dataset(
            {
                'a': (['x'], np.random.rand(4))
            },
            coords = {
                'x': np.arange(4),
            })

    infiles = split_file(tmpdir, d, {'x': 2})

    outpath = tmpdir.join('out.nc')
    c = run_collate(infiles, outpath)

    assert (c.a.data == d.a.data).all()
    assert (c.x.data == d.x.data).all()


def test_split_off_boundary(tmpdir):
    d = xarray.Dataset(
            {
                'a': (['x'], np.random.rand(6))
            },
            coords = {
                'x': np.arange(6),
            })

    infiles = split_file(tmpdir, d, {'x': 3})

    outpath = tmpdir.join('out.nc')
    c = run_collate(infiles, outpath)

    assert (c.a.data == d.a.data).all()
    assert (c.x.data == d.x.data).all()
    
def test_different_compression(tmpdir):
    d = xarray.Dataset(
            {
                'a': (['x'], np.random.rand(4))
            },
            coords = {
                'x': np.arange(4),
            })

    split = split_dataset(d, {'x': 2})

    infiles = [str(x) for x in [tmpdir.join('000.nc'), tmpdir.join('001.nc')]]

    # Save with different compression settings
    split[0].to_netcdf(infiles[0],
            encoding = {'a': {
                'chunksizes': (2,),
                'zlib': True,
                'shuffle': True,
                'complevel': 4,
                }})
    split[1].to_netcdf(infiles[1],
            encoding = {'a': {
                'chunksizes': (2,),
                'zlib': True,
                'shuffle': False,
                'complevel': 6,
                }})

    outpath = tmpdir.join('out.nc')

    c = run_collate(infiles, outpath)

    np.testing.assert_array_equal(d.a, c.a)
    np.testing.assert_array_equal(d.x, c.x)


def test_compression_override(tmpdir):
    d = xarray.Dataset(
            {
                'a': (['x'], np.random.rand(4))
            },
            coords = {
                'x': np.arange(4),
            })

    infiles = split_file(tmpdir, d, {'x': 2})

    outpath = tmpdir.join('out.nc')
    c = run_collate(infiles, outpath, args=['-d','8','--no-shuffle'])

    assert c.a.encoding['complevel'] == 8
    assert c.a.encoding['shuffle'] == False
    np.testing.assert_array_equal(d.a, c.a)

def test_clobber(tmpdir):
    d = xarray.Dataset(
            {
                'a': (['x'], np.random.rand(4))
            },
            coords = {
                'x': np.arange(4),
            })

    infiles = split_file(tmpdir, d, {'x': 2})

    outpath = tmpdir.join('out.nc')
    c = run_collate(infiles, outpath)
    c.close()

    with pytest.raises(subprocess.CalledProcessError):
        c = run_collate(infiles, outpath)

    c = run_collate(infiles, outpath, args=['-f'])

def test_clean(tmpdir):
    d = xarray.Dataset(
            {
                'a': (['x'], np.random.rand(4))
            },
            coords = {
                'x': np.arange(4),
            })

    infiles = split_file(tmpdir, d, {'x': 2})

    outpath = tmpdir.join('out.nc')
    c = run_collate(infiles, outpath, args=['-r'])

    assert not os.path.isfile(infiles[0])
    assert not os.path.isfile(infiles[1])
    assert os.path.isfile(outpath)

def test_1degree(tmpdir):

    # This test dataset is masked, has a missing tile, as there is
    # no data in that tile, and has inconsistent chunk sizes
    testdir = 'test_data'
    testdatadir = os.path.join(testdir, 'output_1deg_masked')

    infiles = glob.glob(os.path.join(testdatadir,'ocean_month.nc.????'))

    # Open file collated by mppnccombine
    d = xarray.open_dataset(os.path.join(testdatadir,'ocean_month.nc'),decode_times=False)

    # No compression
    outpath = tmpdir.join('out.nc')
    c = run_collate(infiles, outpath)

    assert d.equals(c)

    # Compression + shuffle
    outpath = tmpdir.join('out_shuff.nc')
    c = run_collate(infiles, outpath, args=['--shuffle','--deflate','5'])

    assert d.equals(c)

def test_1degree_nc3(tmpdir):

    # This test dataset is masked, has a missing tile, as there is
    # no data in that tile, and has inconsistent chunk sizes
    testdir = 'test_data'
    testdatadir = os.path.join(testdir, 'output_1deg_masked')

    infiles = glob.glob(os.path.join(testdatadir,'ocean_month.nc.????'))

    # Open file collated by mppnccombine
    d = xarray.open_dataset(os.path.join(testdatadir,'ocean_month.nc'),decode_times=False)

    # Convert files to netCDF Classic format and test
    testdir = os.path.join(tmpdir,testdatadir+"_nc3")
    os.makedirs(testdir)
    run_nccopy(['-3'],infiles,testdir)

    infiles = glob.glob(os.path.join(testdir,'ocean_month.nc.????'))

    # No compression
    outpath = tmpdir.join('out.nc')
    c = run_collate(infiles, outpath)

    assert d.equals(c)

    # Compression + shuffle
    outpath = tmpdir.join('out_shuff.nc')
    c = run_collate(infiles, outpath, args=['--shuffle','--deflate','5'])

    assert d.equals(c)

def test_many_files(tmpdir):
    d = xarray.Dataset(
            {
                'a': (['x'], np.random.rand(4000))
            },
            coords = {
                'x': np.arange(4000),
            })

    infiles = split_file(tmpdir, d, {'x': 2})

    outpath = tmpdir.join('out.nc')
    c = run_collate([tmpdir.join('*.nc')], outpath)

    assert d.equals(c)
