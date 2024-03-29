#   You can redistribute this file and/or modify it under the terms of the GNU
#   General Public License GPLv3 (or later), as published by the Free Software
#   Foundation. This file is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.

#   Authors: M. Carla'

import numpy as np
import pyopencl as cl
import pyopencl.tools
import pyopencl.array
from .constants import *
from . import methods
from . import cl_utils
from .methods import FIVED, DOUBLE_PRECISION

import ufo

class Track():
    """
    Track a bunch of particles trhough a ring, allowing to sample particles
    coordinates at different positions at every turn

    line       : Line    -> Ring to be used for tracking
    flags      : int     -> Combination of: LINEAR, FIVED, EXACT, KICK,
                            RADIATION, DOUBLE_PRECISION, ACHROMATIC
    turns      : int     -> Maximum number of turns to be tracked
    particles  : int     -> Number of particles in the bunch
    parameters : list    -> A list of parameters that will not be hardcoded in
                            the pass method. Parameter buffer is accessed via
                            the attribute 'parameters' and should be
                            initialized before tracking
    where      : list    -> List of postions where particles coordinates will
                            be sampled
    dp         : float   -> Bunch relative energy deviation, required only when
                            the flag 'FIVED' is set
    context    : Context -> OpenCL context, as returned by ufo.context() 
    options    : str     -> OpenCL compiler options
    """
    def __init__(self, line, flags=0x00, turns=1000, particles=1000, parameters=[],
                 where=[], dp=0., context=None, options=None):

        self.turns = turns
        self.ctx   = context if context else ufo.DEFAULT_CONTEXT
        self.queue = cl.CommandQueue(self.ctx)
        options    = options if options else ufo.DEFAULT_CL_OPTIONS
        fived      = flags & FIVED
        where      = where if type(where) is list else [where]

        if flags & DOUBLE_PRECISION:
            ufloat = np.float64
        else:
            options += " -cl-single-precision-constant "
            ufloat = np.float32

        self.pool_idx       = np.empty(1, np.int32) #Keep track of already simulated particles
        self.parameters     = np.zeros([particles, len(parameters)], dtype=ufloat)
 
        dump  = "tracks[idx][track_idx][0] = x;\n"
        dump += "tracks[idx][track_idx][1] = px;\n"
        dump += "tracks[idx][track_idx][2] = y;\n"
        dump += "tracks[idx][track_idx][3] = py;\n"
        if not fived:
            dump += "tracks[idx][track_idx][4] = z;\n"
            dump += "tracks[idx][track_idx][5] = dp;\n"

        dump += "track_idx++;\n\n"
 
        methods = line.method(fracture=where, flags=flags, parameters=parameters)
        code = dump.join(methods) #Interleave 'dump' between methods
        self.count = len(where)
        self.where = where #WARNING: dopo lo riassegno, forse posso toglierlo

        if -1 in where:
            code += dump

        flat_line = line.flatten()
        self.where = [(idx, flat_line[int(idx)]) for idx in where] #Warning this is not so nice!!!

        count = self.count * turns
        if count < 1:
            raise Exception("No valid observation points found in 'where'")

        kernel  = f"#define ufloat {'double' if flags & DOUBLE_PRECISION else 'float'}\n"
        kernel += f"#define PARAMETERS_COUNT {len(parameters)}\n"

        kernel += f"__kernel void run(__global ufloat pool[][PARAMETERS_COUNT], __global int *pool_idx, __global ufloat tracks[][{count}][6]) {{\n"
        kernel +=  "    __global ufloat *parameters;\n"
        kernel +=  "    ufloat x, y, z, px, py, dp;\n"
        kernel +=  "    ufloat oodppo;\n" #One Over DP Plus One -> 1/(dp+1)
        kernel +=  "    int idx = get_global_id(0);\n" #First particle has id = thread id
        kernel +=  "    int track_idx;\n"
        kernel +=  "    while (1) {\n" #Loop over particles
        kernel +=  "        parameters = pool[idx];\n"

        for p in ['x', 'px', 'y', 'py', 'z', 'dp']:
            aux = f"parameters[{parameters.index(p)}]" if p in parameters else 0
            kernel += f"        {p}  = {aux};\n"
        if fived:
            kernel += f"        z  = 0.;\n"
            kernel += f"        dp = {dp};\n"

        kernel +=  "        oodppo = 1. / (dp + 1.);\n"
        kernel +=  "        track_idx = 0;\n"
        kernel += f"        for (int turn = 0; turn < {self.turns}; turn++) {{\n"
        kernel += code
        kernel +=  "            if (isnan(x))\n" #Check for lost particles
        kernel +=  "                break;\n"
        kernel +=  "        }\n"
        kernel +=  "        idx = atomic_inc(pool_idx);\n" #Get the next particle
        kernel += f"        if (idx >= {self.parameters.shape[0]})\n"
        kernel +=  "            break;\n" #No left particles
        kernel +=  "    }\n"
        kernel += "}\n"

        self.tracks = np.zeros([particles, count, 6], dtype=ufloat)
        self.src = kernel
        self.kernel = cl.Program(self.ctx, kernel).build(options=options)

        mf = cl.mem_flags
        self.__dev_pool_idx = cl.Buffer(self.ctx, mf.READ_WRITE | mf.COPY_HOST_PTR, hostbuf=self.pool_idx)
        self.__dev_pool     = cl.Buffer(self.ctx, mf.READ_ONLY  | mf.COPY_HOST_PTR, hostbuf=self.parameters)
        self.__dev_tracks   = cl.Buffer(self.ctx, mf.WRITE_ONLY, self.tracks.nbytes)

    def run(self, threads=1):
        """
        Submit the job and start the tracking

        threads : int -> Number of threads to be used for tracking. A number
                         higher than 'particles' may cause undefined behaviour
        """
        self.pool_idx[0] = threads

        cl.enqueue_copy(self.queue, self.__dev_pool_idx, self.pool_idx) #Copy to device
        cl.enqueue_copy(self.queue, self.__dev_pool,     self.parameters)

        evt = self.kernel.run(self.queue, (threads, ), None, self.__dev_pool, self.__dev_pool_idx, self.__dev_tracks) #Run the kernel

        cl.enqueue_copy(self.queue, self.tracks, self.__dev_tracks) #Copy back results

