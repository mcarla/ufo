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

amoeba = '''
#define ALPHA 1.0
#define BETA 0.5
#define GAMMA 2.0

#define get_idx(i, j) (j + i * 4)

void
simplex_sum(ufloat *simplex, ufloat *sum)
{
    int i, j;

    for (j = 0; j < 4; j++)
        for (i = 0, sum[j] = 0.0; i < 4 + 1; i++)
            sum[j] += simplex[get_idx(i, j)];
}

ufloat
amoeba_try(ufloat *simplex, ufloat *y, ufloat *sum, int ihi, ufloat fac, __global ufloat *parameters)
{
    int j;
    ufloat fac1, fac2, y_try;
    ufloat try[4];

    fac1 = (1.0 - fac) / 4;
    fac2 = fac1 - fac;

    for (j = 0; j < 4; j++)
        try[j] = sum[j] * fac1 - simplex[get_idx(ihi, j)] * fac2;

    y_try = method(try, parameters);
    if (y_try < y[ihi]) {
        y[ihi] = y_try;
        for (j = 0; j < 4; j++) {
            sum[j] += try[j] - simplex[get_idx(ihi, j)];
            simplex[get_idx(ihi, j)] = try[j];
        }
    }
    return y_try;
}

int
amoeba(ufloat *simplex, __global ufloat *parameters)
{
    int i, j;
    int ilo, ihi, inhi; /* Low hi and next-hi points */
    int iter;

    ufloat y_try, y_save;
    ufloat tmp_div;
    ufloat y[4 + 1];
    ufloat psum[4];

    iter = 0;
    for (i = 0; i < 4 + 1; i++)
        y[i] = method(&simplex[get_idx(i, 0)], parameters);

    simplex_sum(simplex, psum);
    while (1) { /* Main loop */
        ilo = 1; /* Determine ilo, ihi and inhi */
        ihi = y[1] > y[2] ? (inhi = 2, 1) : (inhi = 1, 2);
        for (i = 0; i < 4 + 1; i++) {
            if (y[i] < y[ilo])
                ilo = i;
            if (y[i] > y[ihi]) {
                inhi = ihi;
                ihi = i;
            } else if
                (y[i] > y[inhi])
                if (i != ihi)
                    inhi = i;
        }
        tmp_div = fabs(y[ihi]) + fabs(y[ilo]);
        if (tmp_div < 1.0e-11)
            tmp_div = 1.0e-11; /* Prevent dividing by zero */

        if (iter >= ITERATIONS)
            break; /* Max iteraions reached */

        y_try = amoeba_try(simplex, y, psum, ihi, -ALPHA, parameters);
        iter++;
        if (y_try <= y[ilo])
            y_try = amoeba_try(simplex, y, psum, ihi, GAMMA, parameters);

        else if (y_try >= y[inhi]) {
            y_save = y[ihi];
            y_try = amoeba_try(simplex, y, psum, ihi, BETA, parameters);
            if (y_try >= y_save) {
                for (i = 0; i < 4 + 1; i++) {
                    if (i != ilo) {
                        for (j = 0; j < 4; j++) {
                            psum[j] = simplex[get_idx(i,   j)] +
                                      simplex[get_idx(ilo, j)];
                            psum[j] *= 0.5;
                            simplex[get_idx(i, j)] = psum[j];
                        }
                        y[i] = method(psum, parameters);
                    }
                }
                iter += 4;
                simplex_sum(simplex, psum);
            }
        }
    }
    return ilo;
}
'''

class ClosedOrbit():
    """
    Find the closed orbit by minimizing the residuals between initial and final
    coordinate of a particle tracked over one turn. The minimization is carried
    out using the Nelder–Mead method.

    line       : Line    -> Ring to be used for tracking
    flags      : int     -> Combination of: LINEAR, FIVED, EXACT, KICK,
                            RADIATION, DOUBLE_PRECISION, ACHROMATIC
    particles  : int     -> Number of particles in the bunch
    parameters : list    -> A list of parameters that will not be hardcoded in
                            the pass method. Parameter buffer is accessed via
                            the attribute 'parameters' and should be
                            initialized before tracking
    dp         : float   -> Bunch relative energy deviation, required only when
                            the flag 'FIVED' is set
    iterations : int     -> Maximum number of minimizer iterations
    context    : Context -> OpenCL context, as returned by ufo.context() 
    options    : str     -> OpenCL compiler options
    """
    def __init__(self, line, flags=0x00, particles=1000, parameters=[], dp=0.,
                 iterations=200, context=cl_utils.context(0), options=None):

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
        self.orbits         = np.zeros([particles, 4], dtype=ufloat)

        mf = cl.mem_flags
        self.__dev_pool_idx = cl.Buffer(self.ctx, mf.READ_WRITE | mf.COPY_HOST_PTR, hostbuf=self.pool_idx)
        self.__dev_pool     = cl.Buffer(self.ctx, mf.READ_ONLY  | mf.COPY_HOST_PTR, hostbuf=self.parameters)
        self.__dev_orbits   = cl.Buffer(self.ctx, mf.WRITE_ONLY, self.orbits.nbytes)

        method = ''.join(line.method(flags=flags, parameters=parameters))

        kernel  = f"#define ufloat {'double' if flags & DOUBLE_PRECISION else 'float'}\n"
        kernel += f"#define PARAMETERS_COUNT {len(parameters)}\n"
        kernel += f"#define ITERATIONS {iterations}\n"

        kernel += 'ufloat method(ufloat orbit[4], __global ufloat *parameters) {'
        kernel += '    ufloat residuals;'
        kernel += '    ufloat x  = orbit[0];'
        kernel += '    ufloat px = orbit[1];'
        kernel += '    ufloat y  = orbit[2];'
        kernel += '    ufloat py = orbit[3];'

        kernel += '    ufloat z  = 0.;'
        kernel += '    ufloat dp = 0.;'
        kernel += '    ufloat oodppo = 1.;'

        kernel += method

        kernel += '    residuals  = (orbit[0] - x ) * (orbit[0] - x );'
        kernel += '    residuals += (orbit[1] - px) * (orbit[1] - px);'
        kernel += '    residuals += (orbit[2] - y ) * (orbit[2] - y );'
        kernel += '    residuals += (orbit[3] - py) * (orbit[3] - py);'
        kernel += '    return residuals;'
        kernel += '}'

        kernel += amoeba

        kernel += f"__kernel void run(__global ufloat pool[][PARAMETERS_COUNT], __global int *pool_idx, __global ufloat orbits[][4]) {{\n"
        kernel +=  "    __global ufloat *parameters;\n"
        kernel +=  "    ufloat simplex[5][4];\n"
        kernel +=  "    ufloat orbit[4];\n"
        kernel +=  "    ufloat oodppo;\n" #One Over DP Plus One -> 1/(dp+1)
        kernel +=  "    int idx = get_global_id(0);\n" #First particle has id = thread id
        kernel +=  "    while (1) {\n" #Loop over particles
        kernel +=  "        parameters = pool[idx];\n"

        #Initialize simplex, if coordinate are provided they will be used as a vertex of the simplex
        for i, p in enumerate(['x', 'px', 'y', 'py']):
            aux = f"parameters[{parameters.index(p)}]" if p in parameters else 0.
            kernel += f"        simplex[0][{i}] = {aux};\n"
            kernel += f"        simplex[1][{i}] = {aux};\n"
            kernel += f"        simplex[2][{i}] = {aux};\n"
            kernel += f"        simplex[3][{i}] = {aux};\n"
            kernel += f"        simplex[4][{i}] = {aux};\n"

        kernel += f"        simplex[1][0] += 1e-4;\n" #Generate the remaining
        kernel += f"        simplex[2][1] += 1e-4;\n" #vertices, by adding a
        kernel += f"        simplex[3][2] += 1e-4;\n" #step in the 4 directions
        kernel += f"        simplex[4][3] += 1e-4;\n"

        kernel +=  "        int ilo = amoeba(simplex, parameters);\n" #Get the next particle

        kernel +=  "        orbits[idx][0] = simplex[ilo][0];\n"
        kernel +=  "        orbits[idx][1] = simplex[ilo][1];\n"
        kernel +=  "        orbits[idx][2] = simplex[ilo][2];\n"
        kernel +=  "        orbits[idx][3] = simplex[ilo][3];\n"

        kernel +=  "        idx = atomic_inc(pool_idx);\n" #Get the next particle
        kernel += f"        if (idx >= {self.parameters.shape[0]})\n"
        kernel +=  "            break;\n" #No left particles
        kernel +=  "    }\n"
        kernel += "}\n"

        self.src = kernel
        self.kernel = cl.Program(context, kernel).build(options=options)

    def run(self, threads=1):
        """
        Submit the job and start closed orbit computation

        threads : int -> Number of threads to be used for the computation. A number
                         higher than 'particles' may cause undefined behaviour
        """
        self.pool_idx[0] = threads

        cl.enqueue_copy(self.queue, self.__dev_pool_idx, self.pool_idx) #Copy to device
        cl.enqueue_copy(self.queue, self.__dev_pool,     self.parameters)

        evt = self.kernel.run(self.queue, (threads, ), None, self.__dev_pool, self.__dev_pool_idx, self.__dev_orbits) #Run the kernel

        cl.enqueue_copy(self.queue, self.orbits, self.__dev_orbits) #Copy back results

