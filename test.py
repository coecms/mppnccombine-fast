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

# Run the collation program
def run_collate(inputs, output, np=2, args=[]):
    try:
        s = subprocess.check_output(
                ['mpirun', '-n', '%d'%np, './mppnccombine-fast', '-o', str(output)] + inputs + args,
                stderr=subprocess.STDOUT)
        print(s.decode('utf-8'))
    except subprocess.CalledProcessError as e:
        print(e.stdout.decode('utf-8'))
        raise
    return xarray.open_dataset(str(output), engine='netcdf4')

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
        infiles.append(str(tmpdir.join('%03d.nc'%i)))

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

    # Fails because of different compressions
    with pytest.raises(subprocess.CalledProcessError) as errinfo:
        c = run_collate(infiles, outpath)

    assert "attributes don't match" in errinfo.value.output.decode()

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
    c = run_collate(infiles, outpath, args=['-d','8','-s'])

    assert c.a.encoding['complevel'] == 8
    assert c.a.encoding['shuffle'] == False

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
