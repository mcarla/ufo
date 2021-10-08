''' UFO is a GPU oriented particle accelerator tracking code '''

from ufo.lattice import *
from ufo.elements import *
from ufo.constants import *
from ufo.track import Track
from ufo.optics import Optics
from ufo.sa import StableAperture
from ufo.co import ClosedOrbit
from ufo.radiation import emittance, RadiationIntegrals
from ufo.chroma import chromaticity
from ufo.rdt import rdt
from ufo.cl_utils import list_devices, context
from ufo.methods import LINEAR, FIVED, EXACT, KICK, RADIATION, DOUBLE_PRECISION, ACHROMATIC

DEFAULT_QUADRUPOLE_SLICES = 8
DEFAULT_BEND_SLICES       = 8
DEFAULT_MULTIPOLE_SLICES  = 4
DEFAULT_SEXTUPOLE_SLICES  = 4
DEFAULT_OCTUPOLE_SLICES   = 4

DEFAULT_CL_OPTIONS = '-cl-fast-relaxed-math -cl-mad-enable'
DEFAULT_CONTEXT = ufo.context(0)

__all__ = ['Lattice', 'Beam', 'Line', 'Track']

try:
    __IPYTHON__
    print(r"""

                            U.F.O.
  

                          \__/  ^__^
                          (oo)  (oo)
                         _/--\  /--\_
                   _.--===0=0====0=0===--._
                  (________________________)
                       /  \________/  \
                      /                \


       An Unreliable, but (Undoubtedly) Fast Optics code
                         Version 0.00

       """)

except:
    pass

