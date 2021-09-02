#   You can redistribute this file and/or modify it under the terms of the GNU
#   General Public License GPLv3 (or later), as published by the Free Software
#   Foundation. This file is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.

#   Authors: M. Carla'

from ufo.lattice import *
from ufo.constants import *
from ufo.track import Track
from ufo.optics import Optics
from ufo.sa import StableAperture
from ufo.radiation import emittance, RadiationIntegrals
from ufo.chroma import chromaticity
from ufo.rdt import rdt
from ufo.cl_utils import list_devices, context
from ufo.methods import LINEAR, FIVED, EXACT, KICK, RADIATION, DOUBLE_PRECISION, ACHROMATIC

DEFAULT_QUADRUPOLE_SLICES = 8
DEFAULT_BEND_SLICES       = 8
DEFAULT_MULTIPOLE_SLICES  = 4
DEFAULT_SEXTUPOLE_SLICES  = 4

DEFAULT_CL_OPTIONS = '-cl-fast-relaxed-math -cl-mad-enable'
DEFAULT_CONTEXT = context(0)
