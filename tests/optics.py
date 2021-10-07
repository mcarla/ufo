import ufo
import color

def test():
    #MAD-X Twiss command results:
    qx = 0.15543311
    qy = 0.362008337 

    optics = ufo.Lattice(path='../optics/alba.mad')

    optics = ufo.Optics(optics.RING, where=[])

    optics.run(threads=1)

    delta = optics.qx - qx
    c = color.OK if abs(delta) < 1e-5 else color.FAIL    
    print(f'Qx -> UFO: {optics.qx:+1.8e}, MAD-X: {qx:+1.8e} ---> Delta: {c} {delta:+1.8e} {color.END}')

    delta = optics.qy - qy
    c = color.OK if abs(delta) < 1e-5 else color.FAIL    
    print(f'Qy -> UFO: {optics.qy:+1.8e}, MAD-X: {qy:+1.8e} ---> Delta: {c} {delta:+1.8e} {color.END}')

