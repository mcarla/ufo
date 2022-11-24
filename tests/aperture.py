import numpy as np
import ufo
import color

def test():
    line = ufo.Line('a_line')
    line.append(ufo.Aperture('ap', window='x > 0.1'))

    parameters=['x']
    tr = ufo.Track(line, turns=1, particles=2, where=[-1], flags=ufo.FIVED, parameters=parameters)

    tr.parameters[0, 0] =  0.0
    tr.parameters[1, 0] =  0.2

    tr.run(threads=1)

    c = color.OK if tr.tracks[0, 0, 0] == 0.0 else color.FAIL
    print(f"{c}Particle 1 should not be lost {color.END}")

    c = color.OK if np.isnan(tr.tracks[1, 0, 0]) else color.FAIL
    print(f"{c}Particle 2 should be lost {color.END}")

