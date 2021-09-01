#   You can redistribute this file and/or modify it under the terms of the GNU
#   General Public License GPLv3 (or later), as published by the Free Software
#   Foundation. This file is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.

#   Authors: M. Carla'

import numpy as np
from . import lattice

def chromaticity(optics):
    dqx0 = dqy0 = dqx = dqy = 0.

    for (_, element), ax, ay, bx, by, dx, dpx in zip(optics.where, optics.ax, optics.ay, optics.bx, optics.by, optics.dx, optics.dpx):
        l = getattr(element, 'length', 0.)
        k = getattr(element, 'k1', 0.)

        if l != 0 and k != 0: #Contribution from quadrupoles
            k2 = np.sqrt(np.absolute(k))
            k2l = l * k2

            gx = (1. + ax*ax) / bx
            gy = (1. + ay*ay) / by

            S  = np.sin (2.*k2l) / 4. / k2
            SH = np.sinh(2.*k2l) / 4. / k2
            C  = (np.cos (k2l)**2 - 1.) / k2 / k2
            CH = (np.cosh(k2l)**2 - 1.) / k2 / k2

            if k > 0:
                dqx0 -= (bx*(l/2. + S ) + ax*C  + gx*(l/2. - S )/k2/k2) * k
                dqy0 += (by*(l/2. + SH) - ay*CH - gy*(l/2. - SH)/k2/k2) * k
            else:
                dqy0 += (by*(l/2. + S ) + ay*C  + gy*(l/2. - S )/k2/k2) * k
                dqx0 -= (bx*(l/2. + SH) - ax*CH - gx*(l/2. - SH)/k2/k2) * k

        k = getattr(element, 'k2', 0.)
        if l != 0 and k != 0: #Contribution from sextupoles
            gx = (1. + ax*ax) / bx
            gy = (1. + ay*ay) / by

            dqx += (dx * (bx + (gx*l/3. - ax)*l) + dpx*l*(bx/2. - l*(2.*ax/3. + gx*l/4.))) * l * k
            dqy -= (dx * (by + (gy*l/3. - ay)*l) + dpx*l*(by/2. - l*(2.*ay/3. + gy*l/4.))) * l * k

    dqx0 /= 4. * np.pi
    dqy0 /= 4. * np.pi
    dqx  /= 4. * np.pi
    dqy  /= 4. * np.pi

    dqx += dqx0
    dqy += dqy0

    return (dqx0, dqy0), (dqx, dqy)

