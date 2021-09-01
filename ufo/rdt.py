#   You can redistribute this file and/or modify it under the terms of the GNU
#   General Public License GPLv3 (or later), as published by the Free Software
#   Foundation. This file is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.

#   Authors: M. Carla'

import numpy as np
from . import lattice

def rdt(optics):
    f3000 = f1200 = f1020 = f0120 = f0111 = 0.

    for (_, element), bx, by, mux, muy in zip(optics.where, optics.bx, optics.by, optics.mux, optics.muy):
        l = getattr(element, 'length', 0.)
        k = getattr(element, 'k2', 0.)

        if l != 0 and k != 0:
            klb = k * l * np.sqrt(bx)

            f3000 += klb * bx * np.exp( 3j *   mux)
            f1200 += klb * bx * np.exp(-1j *   mux)
            f1020 += klb * by * np.exp( 1j * ( mux + 2.*muy))
            f0120 += klb * by * np.exp( 1j * (-mux + 2.*muy))
            f0111 += klb * by * np.exp(-1j *   mux)

    return f3000, f1200, f1020, f0120, f0111

