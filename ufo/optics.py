#   You can redistribute this file and/or modify it under the terms of the GNU
#   General Public License GPLv3 (or later), as published by the Free Software
#   Foundation. This file is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.

#   Authors: M. Carla'

import numpy as np
import math
import pyopencl as cl
import pyopencl.tools
from .constants import *
from . import methods
from . import mad
from . import cl_utils
from .methods import FIVED, LINEAR, ACHROMATIC, DOUBLE_PRECISION
from . import track

import ufo

class Optics(track.Track):
    """
    Compute betatron and dispersion function of a given line with periodic or
    user specified boundary conditions.

    line       : Line    -> Ring to be used for tracking
    periodic   : Boolen  -> If True periodic boundary conditions are used.
                            If False user specified initial condition are used,
                            initial conditions should be specified before
                            calling run() by means of the attributes: ax, ay,
                            bx, by, dx, dy, dpx, dpy.

    flags      : int     -> Combination of: LINEAR, FIVED, EXACT, KICK,
                            RADIATION, DOUBLE_PRECISION, ACHROMATIC
    parameters : list    -> A list of parameters that will not be hardcoded in
                            the pass method. Parameter buffer is accessed via
                            the attribute 'parameters' and should be
                            initialized before tracking
    where      : list    -> List of postions where optics functions coordinates
                            will be evaluated
    context    : Context -> OpenCL context, as returned by ufo.context() 
    options    : str     -> OpenCL compiler options
    """
    def __init__(self, line, where=[], periodic=True, parameters=[], flags=LINEAR | ACHROMATIC,
                 context=cl_utils.context(0), options=None):

        self.ctx      = context if context else ufo.DEFAULT_CONTEXT
        self.queue    = cl.CommandQueue(self.ctx)
        self.periodic = periodic
        options       = options if options else ufo.DEFAULT_CL_OPTIONS
        fived         = flags & FIVED
        where         = where if type(where) is list else [where]
        particles     = 3

        if periodic and -1 not in where and line[-1] not in where:
            where.append(-1) #One-turn-map is required to calculate the periodic solutions

        if flags & DOUBLE_PRECISION:
            ufloat = np.float64
        else:
            options += " -cl-single-precision-constant "
            ufloat = np.float32

        self.pool_idx   = np.empty(1, np.int32) #Keep track of already simulated particles
        self.parameters = np.zeros([1, max(1, len(parameters))], dtype=ufloat) #If set to 0 size cl.Buffer would complain!

        dump  = "tracks[particle_idx][track_idx][0] = x;\n"
        dump += "tracks[particle_idx][track_idx][1] = px;\n"
        dump += "tracks[particle_idx][track_idx][2] = y;\n"
        dump += "tracks[particle_idx][track_idx][3] = py;\n"
        dump += "track_idx++;\n\n"

        methods = line.method(fracture=where, flags=flags, parameters=parameters)
        code = dump.join(methods) #Interleave 'dump' between methods
        self.count = len(where) - 1

        if -1 in where:
            self.count += 1
            code += dump

        flat_line = line.flatten()
        self.where = [(idx, flat_line[int(idx)]) for idx in where] #Warning this is not so nice!!!

        if self.count < 1:
            raise Exception("No valid observation points found in 'where'")

        kernel  = f"#define ufloat {'double' if flags & DOUBLE_PRECISION else 'float'}\n"
        kernel += f"__kernel void run(__global ufloat parameters[][{len(parameters)}], __global int *pool_idx, __global ufloat tracks[][{self.count}][6]) {{\n"
        kernel +=  "    ufloat x, y, z, px, py, dp;\n"
        kernel +=  "    ufloat oodppo;\n" #One Over DP Plus One -> 1/(dp+1)
        kernel +=  "    int particle_idx = get_global_id(0);\n" #First particle has id = thread id
        kernel +=  "    int idx = 0;\n"
        kernel +=  "    int track_idx;\n"
        kernel +=  "    while (1) {\n" #Loop over particles

        kernel += f"        x = y = z = px = py = dp = 0.;\n"
        kernel += f"        if (particle_idx == 0) {{x  = 1.; y  = 1.;}}\n"
        kernel += f"        if (particle_idx == 1) {{px = 1.; py = 1.;}}\n"
        kernel += f"        if (particle_idx == 2) dp = 1.;\n"

        if 'dp' in parameters:
            kernel += f"        dp += parameters[0][{parameters.index('dp')}];\n"

        kernel +=  "        oodppo = 1. / (dp + 1.);\n"
        kernel +=  "        track_idx = 0;\n"
        kernel += code
        kernel +=  "        particle_idx = atomic_inc(pool_idx);\n" #Get the next particle
        kernel += f"        if (particle_idx >= {particles})\n"
        kernel +=  "            break;\n" #No left particles
        kernel +=  "    }\n"
        kernel += "}\n"

        self.track = np.zeros([particles, self.count, 6], dtype=ufloat)
        self.src = kernel
        self.kernel = cl.Program(context, kernel).build(options=options)

        mf = cl.mem_flags
        self.__dev_pool_idx = cl.Buffer(self.ctx, mf.READ_WRITE | mf.COPY_HOST_PTR, hostbuf=self.pool_idx)
        self.__dev_pool     = cl.Buffer(self.ctx, mf.READ_ONLY  | mf.COPY_HOST_PTR, hostbuf=self.parameters)
        self.__dev_track    = cl.Buffer(self.ctx, mf.WRITE_ONLY, self.track.nbytes)

        self.ax0 = self.bx0 = self.dx0 = self.dpx0 = 0.
        self.ay0 = self.by0 = self.dy0 = self.dpy0 = 0.

        self.ax  = np.empty(self.count)
        self.ay  = np.empty(self.count)
        self.bx  = np.empty(self.count)
        self.by  = np.empty(self.count)
        self.dx  = np.empty(self.count)
        self.dy  = np.empty(self.count)
        self.dpx = np.empty(self.count)
        self.dpy = np.empty(self.count)
        self.mux = np.empty(self.count)
        self.muy = np.empty(self.count)

    def run(self, threads=3):
        """
        Submit the job and start the tracking

        threads : int -> Number of threads to be used for tracking. A number
                         higher than '3*number of parameters' may cause undefined behaviour
        """
        self.pool_idx[0] = threads

        cl.enqueue_copy(self.queue, self.__dev_pool_idx, self.pool_idx) #Copy to device
        cl.enqueue_copy(self.queue, self.__dev_pool,     self.parameters)

        evt = self.kernel.run(self.queue, (threads, ), None, self.__dev_pool, self.__dev_pool_idx, self.__dev_track) #Run the kernel

        cl.enqueue_copy(self.queue, self.track, self.__dev_track) #Copy back results

        #M00 = re[x], M10 = re[px], M01 = im[x], M11 = im[px]
        if self.periodic:
            re = self.track[0][-1]
            im = self.track[1][-1]
            dp = self.track[2][-1]

            #Horizontal plane
            cos_nu = (re[0] + im[1]) * 0.5
            if np.absolute(cos_nu) > 1.:
                self.ax0 = self.bx0 = self.qx = np.nan
                return

            #Compute alpha / beta / gamma
            self.bx0 = np.absolute(im[0] / np.sin(np.arccos(cos_nu)))
            sin_nu = im[0] / self.bx0
            nu = np.arctan2(sin_nu, cos_nu)

            if nu < 0.:
                nu += 2.*np.pi
 
            self.ax0 = (re[0] - im[1]) / (2. * np.sin(nu))
            self.qx = nu / (2. * np.pi)

            #Compute dispersion
            k = 2. - re[0]  - im[1]
            self.dx0  = ((1. - im[1]) * dp[0] + im[0] * dp[1]) / k
            self.dpx0 = ((1. - re[0]) * dp[1] + re[1] * dp[0]) / k

            #Vertical plane
            cos_nu = (re[2] + im[3]) * 0.5
            if np.absolute(cos_nu) > 1.:
                self.ay0 = self.by0 = self.qy = np.nan
                return

            #Compute alpha / beta / gamma
            self.by0 = np.absolute(im[2] / np.sin(np.arccos(cos_nu)))
            sin_nu = im[2] / self.by0
            nu = np.arctan2(sin_nu, cos_nu)

            if nu < 0.:
                nu += 2.*np.pi
 
            self.ay0 = (re[2] - im[3]) / (2. * np.sin(nu))
            self.qy = nu / (2. * np.pi)

        for idx, (re, im, dp) in enumerate(zip(self.track[0], self.track[1], self.track[2])):
           #Horizontal plane
           gamma0 = (1. + self.ax0 * self.ax0) / self.bx0
           self.ax[idx]  = -re[0] * re[1] * self.bx0
           self.ax[idx] += (re[0] * im[1] + im[0] * re[1]) * self.ax0
           self.ax[idx] += -im[0] * im[1] * gamma0

           self.bx[idx]  =  re[0] * re[0] * self.bx0
           self.bx[idx] += -re[0] * im[0] * self.ax0 * 2.
           self.bx[idx] +=  im[0] * im[0] * gamma0
           #if self.bx[idx] <= 0.:
           #    raise ArithmeticError('Reached numerical precision limit') 

           self.dx[idx]  =  re[0] * self.dx0 + im[0] * self.dpx0 + dp[0]
           self.dpx[idx] =  re[1] * self.dx0 + im[1] * self.dpx0 + dp[1]

           sin_mu = im[0] / np.sqrt(self.bx0 * self.bx[idx])
           cos_mu = re[0] * np.sqrt(self.bx0 / self.bx[idx]) - self.ax0 * sin_mu
           self.mux[idx] = np.arctan2(sin_mu, cos_mu)

           #Vertical plane
           gamma0 = (1. + self.ay0 * self.ay0) / self.by0
           self.ay[idx]  = -re[2] * re[3] * self.by0
           self.ay[idx] += (re[2] * im[3] + im[2] * re[3]) * self.ay0
           self.ay[idx] += -im[2] * im[3] * gamma0

           self.by[idx]  =  re[2] * re[2] * self.by0
           self.by[idx] += -re[2] * im[2] * self.ay0 * 2.
           self.by[idx] +=  im[2] * im[2] * gamma0
           #if self.by[idx] <= 0.:
           #    raise ArithmeticError('Reached numerical precision limit') 

           self.dy[idx]  =  re[2] * self.dy0 + im[2] * self.dpy0 + dp[2]
           self.dpy[idx] =  re[3] * self.dy0 + im[3] * self.dpy0 + dp[3]

           sin_mu = im[2] / np.sqrt(self.by0 * self.by[idx])
           cos_mu = re[2] * np.sqrt(self.by0 / self.by[idx]) - self.ay0 * sin_mu
           self.muy[idx] = np.arctan2(sin_mu, cos_mu)

