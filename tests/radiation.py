import ufo
import color

def test():
    I1 =  0.289626312
    I2 =  0.8979648783
    I3 =  0.1283331437
    I4 = -0.3327286365
    I5 =  0.0005327915674
    U0 =  1.
    E0 =  4
    alpha_c = 0.0008877928529  

    optics = ufo.Lattice(path='../optics/alba.mad')

    where = optics.RING.find(lambda e: hasattr(e, 'angle'))
    of = ufo.Optics(optics.RING, where=where)

    of.run(threads=1)

    rad = ufo.emittance(of)

    delta = rad.I1 - I1
    c = color.OK if abs(delta) < 1e-2 else color.FAIL
    print(f"I1 -> UFO: {rad.I1:+1.8e}, MAD-X: {I1:+1.8e} ---> Delta: {c} {delta:+1.8e} {color.END}")

    delta = rad.I2 - I2
    c = color.OK if abs(delta) < 1e-2 else color.FAIL
    print(f"I2 -> UFO: {rad.I2:+1.8e}, MAD-X: {I2:+1.8e} ---> Delta: {c} {delta:+1.8e} {color.END}")

    delta = rad.I4 - I4
    c = color.OK if abs(delta) < 1e-2 else color.FAIL
    print(f"I4 -> UFO: {rad.I4:+1.8e}, MAD-X: {I4:+1.8e} ---> Delta: {c} {delta:+1.8e} {color.END}")

    delta = rad.I5 - I5
    c = color.OK if abs(delta) < 1e-2 else color.FAIL
    print(f"I5 -> UFO: {rad.I5:+1.8e}, MAD-X: {I5:+1.8e} ---> Delta: {c} {delta:+1.8e} {color.END}")

    print(f"U0 [MeV/turn] -> UFO: {rad.U0:2.10e}, MAD-X: {'?'}, Delta: {'?'}")
    print(f"Ex [nm*rad] -> UFO: {rad.Ex:2.10e}, MAD-X: {'?'}, Delta: {'?'}")

    alpha_ufo = rad.I1 / optics.RING.length
    delta = alpha_ufo - alpha_c
    c = color.OK if abs(delta) < 1e-2 else color.FAIL
    print(f"alpha_c -> UFO: {alpha_ufo:+1.8e}, MAD-X: {alpha_c:+1.8e}, Delta: {c} {delta:+1.8e} {color.END}")

