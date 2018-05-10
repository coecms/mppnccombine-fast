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
import matplotlib.pyplot as plt


d = xarray.open_dataset('/g/data/w35/saw562/test0.nc', chunks={'time': 1, 'st_ocean': 1})
temp = d.isel(time=0, st_ocean=0).temp
#d.isel(time=0, st_ocean=0).temp.plot.imshow()
temp.plot.imshow()
#plt.imshow(temp.values)
plt.show()
