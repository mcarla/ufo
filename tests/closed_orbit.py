import numpy as np
import ufo
import color

def test():
    lattice = ufo.Lattice(path='../optics/alba.mad')

    lattice.QD1.dy = 1e-4
    lattice.QD1.dx = 1e-4

    orbit = ufo.ClosedOrbit(lattice.RING, particles=1)
    orbit.run()

    #Check against tracking
    parameters = ['x', 'px', 'y', 'py']
    track = ufo.Track(lattice.RING, parameters=parameters, turns=1, particles=1, where=[-1])

    track.parameters[0, 0] = orbit.orbits[0, 0]
    track.parameters[0, 1] = orbit.orbits[0, 1]
    track.parameters[0, 2] = orbit.orbits[0, 2]
    track.parameters[0, 3] = orbit.orbits[0, 3]

    track.run()

    delta = track.tracks[0, 0, :4] - orbit.orbits[0]

    for idx, label in enumerate(parameters):
        c = color.OK if abs(delta[idx]) < 1e-6 else color.FAIL
        print(f"Closed-orbit vs Tracking: Delta-{label}: {c} {delta[idx]:1.8e} {color.END}")

