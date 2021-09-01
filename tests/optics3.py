import ufo
import color

def test():
    #MAD-X Twiss command results:
    dqx_mad = -40.13060686
    dqy_mad = -27.05066939

    delta_dp = 1e-4

    lattice = ufo.Lattice(path='../optics/alba.mad')

    optics = ufo.Optics(lattice.RING, where=[], parameters=['dp'], flags=ufo.LINEAR)

    optics.parameters[0] = 0.0
    optics.run(threads=1)
    qx0 = optics.qx
    qy0 = optics.qy

    optics.parameters[0] = delta_dp
    optics.run(threads=1)
    dqx = (optics.qx - qx0) / delta_dp
    dqy = (optics.qy - qy0) / delta_dp

    delta = dqx - dqx_mad
    c = color.OK if abs(delta) < 1e-1 else color.FAIL    
    print(f'Qx -> UFO: {dqx:+1.8e}, MAD-X: {dqx_mad:+1.8e} ---> Delta: {c} {delta:+1.8e} {color.END}')

    delta = dqy - dqy_mad
    c = color.OK if abs(delta) < 1e-1 else color.FAIL    
    print(f'Qy -> UFO: {dqy:+1.8e}, MAD-X: {dqy_mad:+1.8e} ---> Delta: {c} {delta:+1.8e} {color.END}')

