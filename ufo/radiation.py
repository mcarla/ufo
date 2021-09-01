#   You can redistribute this file and/or modify it under the terms of the GNU
#   General Public License GPLv3 (or later), as published by the Free Software
#   Foundation. This file is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.

#   Authors: M. Carla'

import numpy as np
from .constants import *
from . import lattice

class RadiationIntegrals():
    def __init__(self, I1, I2, I4, I5):
        self.I1 = I1
        self.I2 = I2
        self.I4 = I4
        self.I5 = I5
        self.Jx = 1. - I4 / I2
        self.Jt = 1. + I4 / I2

    def __get_U0(self, beam=lattice.Beam()):
        return radiation_loss_constant * (beam.energy * 1e-9)**4 * self.I2

    def __get_Ex(self, beam=lattice.Beam()):
        k_emitt = quantum_radiation_constant * beam.gamma * beam.gamma
        return k_emitt * self.I5 / self.I2 / self.Jx

    U0  = property(__get_U0,  None)
    Ex  = property(__get_Ex,  None)

def emittance(optics):
    I1 = I2 = I4 = I5 = 0

    for (_, element), alpha, beta, d, dp in zip(optics.where, optics.ax, optics.bx, optics.dx, optics.dpx):
        # From SLAC 1193 / Elegant
        angle  = getattr(element, 'angle', 0.)
        if angle == 0.:
            continue

        length = getattr(element, 'length', 0.)
        k      = getattr(element, 'k1', 0.)
        k2 = np.sqrt(np.absolute(k))
        k2l = length * k2
        rho = length / angle
        gamma = (1. + alpha * alpha) / beta

        if k < 0: #I should add 1/rho^2 see elegant source
                C = np.cosh(k2l)
                S = np.sinh(k2l)
        else:
                C = np.cos(k2l)
                S = np.sin(k2l)

        I2 += 1. / rho**2. * length

        #Dispersion integral
        Id  = d * S / k2l + dp * (1. - C) / (length * k)
        Id += (k2l - S) / (rho * k2l * k)
        I1 += Id * length / rho
        I4 += Id * (1. / rho**2. + 2.*k) * length / rho

        H  = gamma * d**2. + 2.*alpha * d * dp + dp**2 * beta
        H += 2.*angle * (-(gamma * d + alpha * dp) * (k2l - S) / (k2 * k * length**2)
             + (alpha * d + beta * dp) * (1. - C) / (k * length**2))

        H += angle**2 * (gamma * (3. * k2l - 4. * S + S*C) / (2. * k2 * k**2 * length**3)
             - alpha * (1. - C)**2 / (k**2 * length**3)
             + beta * (k2l - C * S) / (2. * k2l * k * length**2))

        I5 += 1. / np.absolute(rho**3) * H * length

    return RadiationIntegrals(I1, I2, I4, I5)

