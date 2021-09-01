#   You can redistribute this file and/or modify it under the terms of the GNU
#   General Public License GPLv3 (or later), as published by the Free Software
#   Foundation. This file is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.

#   Authors: M. Carla'

import pyopencl as cl
import pyopencl.tools

# Check what OpenCL devices are available
def list_devices(device=None):

    platforms = cl.get_platforms()

    idx = 0
    for platform in cl.get_platforms():
        for d in platform.get_devices():
            if device == idx:
                return d
            elif not device:
                print(f'{idx}:   {d.name}')
            idx += 1
 
def context(device):
    return cl.Context([list_devices(device=device)])

