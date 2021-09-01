import ufo
import color

def test():
    align = ufo.Lattice(path='../optics/align.mad')

    align.QF.dx = 1e-4
    fodo = ufo.Optics(align.FODO, where=[-1])
    equivalent = ufo.Optics(align.EQUIVALENT, where=[-1])

    fodo.run(threads=1)
    equivalent.run(threads=1)

    delta = fodo.dx[0] - equivalent.dx[0]
    c = color.OK if abs(delta) < 1e-6 else color.FAIL
    print(f'Dx -> Align: {fodo.dx[0]:+1.6e}, Equivalent: {equivalent.dx[0]:+1.6e} ---> Delta: {c} {delta:+1.5e} {color.END}')

