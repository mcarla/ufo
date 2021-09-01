import numpy as np
import ufo
import color

def test():
    x = np.array([-0.0008304642239, -0.001747306891, -0.0007418783837, 0.001079715459, 0.001713478007,
                   0.000462185744,  -0.00129757264,  -0.001629827644, -0.0001690544829,0.001477701296])

    parameters=['x']

    optics = ufo.Lattice(path='../optics/fodo.mad')

    tr = ufo.Track(optics.RING, turns=10, particles=2, where=[-1], flags=ufo.FIVED, parameters=parameters)
    tr.parameters[:, 0] =  0.001

    tr.run(threads=1)

    delta = np.max(np.absolute(tr.track[0][:, 0] - x))
    c = color.OK if abs(delta) < 1e-6 else color.FAIL
    print(f"Particle 1: UFO/MAD-X Max |Delta| [m]: {c} {delta:+1.8e} {color.END}")

    delta = np.max(np.absolute(tr.track[1][:, 0] - x))
    c = color.OK if abs(delta) < 1e-6 else color.FAIL
    print(f"Particle 2: UFO/MAD-X Max |Delta| [m]: {c} {delta:+1.8e} {color.END}")

