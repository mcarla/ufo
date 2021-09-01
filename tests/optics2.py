import ufo
import color

def test():
    optics = ufo.Lattice(path='../optics/alba.mad')

    periodic  = ufo.Optics(optics.RING, where=[-1], periodic=True)
    propagate = ufo.Optics(optics.RING, where=[-1], periodic=False)

    periodic.run(threads=1)
    propagate.ax0 = periodic.ax[-1]
    propagate.ay0 = periodic.ay[-1]
    propagate.bx0 = periodic.bx[-1]
    propagate.by0 = periodic.by[-1]
    propagate.dx0 = periodic.dx[-1]
    propagate.dpx0 = periodic.dpx[-1]

    propagate.run(threads=1)

    delta = periodic.ax[-1] - propagate.ax[-1]
    c = color.OK if abs(delta) < 1e-5 else color.FAIL
    print(f'Ax -> Periodic: {periodic.ax[-1]:+1.5e}, Propagate: {propagate.ax[-1]:+1.5e} ---> Delta: {c} {delta:+1.5e} {color.END}')

    delta = periodic.ay[-1] - propagate.ay[-1]
    c = color.OK if abs(delta) < 1e-5 else color.FAIL
    print(f'Ay -> Periodic: {periodic.ay[-1]:+1.5e}, Propagate: {propagate.ay[-1]:+1.5e} ---> Delta: {c} {delta:+1.5e} {color.END}')

    delta = periodic.bx[-1] - propagate.bx[-1]
    c = color.OK if abs(delta) < 1e-5 else color.FAIL
    print(f'Bx -> Periodic: {periodic.bx[-1]:+1.5e}, Propagate: {propagate.bx[-1]:+1.5e} ---> Delta: {c} {delta:+1.5e}, {color.END}')

    delta = periodic.by[-1] - propagate.by[-1]
    c = color.OK if abs(delta) < 1e-5 else color.FAIL
    print(f'By -> Periodic: {periodic.by[-1]:+1.5e}, Propagate: {propagate.by[-1]:+1.5e} ---> Delta: {c} {delta:+1.5e} {color.END}')


    delta = periodic.dx[-1] - propagate.dx[-1]
    c = color.OK if abs(delta) < 1e-5 else color.FAIL
    print(f'Dx  -> Periodic: {periodic.dx[-1]:+1.5e}, Propagate: {propagate.dx[-1]:+1.5e} ---> Delta: {c} {delta:+1.5e} {color.END}')

    delta = periodic.dpx[-1] - propagate.dpx[-1]
    c = color.OK if abs(delta) < 1e-5 else color.FAIL
    print(f'Dpx -> Periodic: {periodic.dpx[-1]:+1.5e}, Propagate: {propagate.dpx[-1]:+1.5e} ---> Delta: {c} {delta:+1.5e} {color.END}')

