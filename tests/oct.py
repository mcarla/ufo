import numpy as np
import ufo
import color


parameters=['x']

optics = ufo.Lattice(path='../optics/oct.mad')

tr = ufo.Track(optics.RING, turns=10, particles=2, where=[-1], flags=ufo.FIVED, parameters=parameters)
tr.parameters[:, 0] =  0.01

tr.run(threads=1)

print(tr.track[0][:, 0])

