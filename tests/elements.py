import ufo
import color

def test():
    elements = ufo.Lattice(path='../optics/elements.mad')
    elements.TEST.method()
    print(f'{color.OK} OK {color.END}')
