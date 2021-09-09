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

class StableAperture():
    """
    Find stable aperture by tracking a bunch of particles, returns the number
    of turns for each particle before the particle getting lost

    line       : Line    -> Ring to be used for tracking
    flags      : int     -> Combination of: LINEAR, FIVED, EXACT, KICK,
                            RADIATION, DOUBLE_PRECISION, ACHROMATIC
    turns      : int     -> Maximum number of turns to be tracked
    particles  : int     -> Number of particles in the bunch
    parameters : list    -> A list of parameters that will not be hardcoded in
                            the pass method. Parameter buffer is accessed via
                            the attribute 'parameters' and should be
                            initialized before tracking
    dp         : float   -> Bunch relative energy deviation, required only when
                            the flag 'FIVED' is set
    context    : Context -> OpenCL context, as returned by ufo.context() 
    options    : str     -> OpenCL compiler options
    """
    def __init__(self, line, flags=0x00, turns=1000, particles=1000, parameters=[], dp=0.,
                 context=cl_utils.context(0), options=None):

        self.turns = turns
        self.ctx   = context if context else ufo.DEFAULT_CONTEXT
        self.queue = cl.CommandQueue(self.ctx)
        options    = options if options else ufo.DEFAULT_CL_OPTIONS
        fived      = flags & FIVED

        if flags & DOUBLE_PRECISION:
            ufloat = np.float64
        else:
            options += " -cl-single-precision-constant "
            ufloat = np.float32

        self.pool_idx       = np.empty(1, np.int32) #Keep track of already simulated particles
        self.parameters     = np.zeros([particles, max(1, len(parameters))], dtype=ufloat) #If set to 0 size cl.Buffer would complain!
        self.lost           = np.empty(particles, dtype=np.uint32) #Keep track of lost particles

        mf = cl.mem_flags
        self.__dev_pool_idx = cl.Buffer(self.ctx, mf.READ_WRITE | mf.COPY_HOST_PTR, hostbuf=self.pool_idx)
        self.__dev_pool     = cl.Buffer(self.ctx, mf.READ_ONLY  | mf.COPY_HOST_PTR, hostbuf=self.parameters)
        self.__dev_lost     = cl.Buffer(self.ctx, mf.WRITE_ONLY, self.lost.nbytes)

        partial_turn = ''
        full_turn    = ''

        '''
        if 's' in parameters: #First partial turn with fractured elements
            count = 0 #Element index
            s_idx = parameters.index('s') #S parameter index
            for element in line.flatten():
                prms = [idx for (idx, p) in enumerate(parameters) if p[0]==element.label]
                prms = {parameters[p][1]: f'parameters[idx][{p}]' for p in prms}
                prms['flags'] = flags

                partial_turn +=  "fracture = 0.;\n" #Fracture point
                partial_turn += f"if ({count} == floor(parameters[idx][{s_idx}])) {{\n" #Check if s falls into the current element

                for p in ['x', 'px', 'y', 'py', 'z', 'dp']: #Reset particle coordinates to initial value
                    aux = f"parameters[idx][{parameters.index(p)}]" if p in parameters else 0
                    partial_turn += f"        {p}  = {aux};\n"
                if fived:
                    partial_turn += f"        z  = 0.;\n"
                    partial_turn += f"        dp = {dp};\n"

                partial_turn += f"    fracture = parameters[idx][{s_idx}] - floor(parameters[idx][{s_idx}]);\n"
                partial_turn +=  "}\n"

                partial_turn += element.method(**prms, fracture='fracture')
                count += 1
        '''
        full_turn = ''.join(line.method(flags=flags, parameters=parameters))

        kernel  = f"#define ufloat {'double' if flags & DOUBLE_PRECISION else 'float'}\n"
        kernel += f"__kernel void run(__global ufloat parameters[][{len(parameters)}], __global int *pool_idx, __global uint *lost) {{\n"
        kernel +=  "    ufloat x, y, z, px, py, dp;\n"
        kernel +=  "    ufloat fracture;\n"
        kernel +=  "    ufloat oodppo;\n" #One Over DP Plus One -> 1/(dp+1)
        kernel +=  "    int turn, idx = get_global_id(0);\n" #First particle has id = thread id
        kernel +=  "    while (1) {\n" #Loop over particles
        kernel +=  "        lost[idx] = 0;\n" #Will stay 0 if stable

        for p in ['x', 'px', 'y', 'py', 'z', 'dp']:
            aux = f"parameters[idx][{parameters.index(p)}]" if p in parameters else 0
            kernel += f"        {p}  = {aux};\n"
        if fived:
            kernel += f"        z  = 0.;\n"
            kernel += f"        dp = {dp};\n"

        kernel +=  "        oodppo = 1. / (dp + 1.);\n"

        if 's' in parameters: #Run the first partial turn if 's' is specified
            kernel += partial_turn

        kernel += f"        for (turn = 0; turn < {self.turns}; turn++) {{\n" 
        kernel +=               full_turn
        kernel +=  "            if (fabs(x) > 1. || fabs(y) > 1. || isnan(x))\n" #Check for lost particles
        kernel +=  "                break;\n"
        kernel +=  "        }\n"
        kernel +=  "        lost[idx] = turn;\n" #Store last stable turn
        kernel +=  "        idx = atomic_inc(pool_idx);\n" #Get the next particle
        kernel += f"        if (idx >= {self.parameters.shape[0]})\n"
        kernel +=  "            break;\n" #No left particles
        kernel +=  "    }\n"
        kernel += "}\n"

        self.src = kernel
        self.kernel = cl.Program(context, kernel).build(options=options)

    def run(self, threads=1):
        """
        Submit the job and start the tracking

        threads : int -> Number of threads to be used for tracking. A number
                         higher than 'particles' may cause undefined behaviour
        """
        self.pool_idx[0] = threads

        cl.enqueue_copy(self.queue, self.__dev_pool_idx, self.pool_idx) #Copy to device
        cl.enqueue_copy(self.queue, self.__dev_pool,     self.parameters)

        evt = self.kernel.run(self.queue, (threads, ), None, self.__dev_pool, self.__dev_pool_idx, self.__dev_lost) #Run the kernel

        cl.enqueue_copy(self.queue, self.lost, self.__dev_lost) #Copy back results

