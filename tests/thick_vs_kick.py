import ufo
import color

def test():
    ufo.DEFAULT_QUADRUPOLE_SLICES = 16
    ufo.DEFAULT_BEND_SLICES       = 48

    optics = ufo.Lattice(path='../optics/alba.mad')

    kick  = ufo.Optics(optics.RING, flags=ufo.LINEAR | ufo.ACHROMATIC | ufo.KICK)
    thick = ufo.Optics(optics.RING, flags=ufo.LINEAR | ufo.ACHROMATIC)

    kick.run(threads=1)
    thick.run(threads=1)

    delta = thick.qx - kick.qx
    c = color.OK if abs(delta) < 5e-4 else color.FAIL
    print(f'Qx -> Thick {thick.qx:1.8f}, Kick: {kick.qx:1.8f} ---> Delta: {c} {delta:+1.6e} {color.END}')

    delta = thick.qy - kick.qy
    c = color.OK if abs(delta) < 5e-4 else color.FAIL
    print(f'Qy -> Thick {thick.qy:1.8f}, Kick: {kick.qy:1.8f} ---> Delta: {c} {delta:+1.6e} {color.END}')

