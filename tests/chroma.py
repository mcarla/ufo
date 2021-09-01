import numpy as np
import ufo
import color

def test():
    #MAD-X Twiss command results:
    dqx = 1.752728636
    dqy = 3.589260652 

    optics = ufo.Lattice(path='../optics/alba.mad')

    where = optics.RING.find(lambda e: hasattr(e, 'k1') or hasattr(e, 'k2'))

    of = ufo.Optics(optics.RING, where=where)

    of.run(threads=1)
    dq = ufo.chromaticity(of)

    delta = dq[1][0] - dqx
    c = color.OK if abs(delta) < 1e-3 else color.FAIL
    print(f'dQx -> UFO: {dq[1][0]:+1.8e}, MAD-X: {dqx:+1.8e} ---> Delta: {c} {delta:+1.8e} {color.END}')


    delta = dq[1][1] - dqy
    c = color.OK if abs(delta) < 1e-3 else color.FAIL
    print(f'dQy -> UFO: {dq[1][1]:+1.8e}, MAD-X: {dqy:+1.8e} ---> Delta: {c} {delta:+1.8e} {color.END}')

