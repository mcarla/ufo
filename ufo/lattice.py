#   You can redistribute this file and/or modify it under the terms of the GNU
#   General Public License GPLv3 (or later), as published by the Free Software
#   Foundation. This file is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.

#   Authors: M. Carla'

import numpy as np
import pyopencl as cl
import pyopencl.tools
from .constants import *
from . import methods
from . import mad
from . import cl_utils
from .methods import FIVED

import ufo

class Element():
    def __init__(self):
        pass

    def __str__(self):
        return self.label 

    def __repr__(self):  
        msg = f'{self.label: <8}: {type(self).__name__.upper(): <12} '
        return msg + '   '.join([f'{p}={self.__dict__[p]}' for p in self.parameters])

    def __replace__(self, fragment, **kwargs):
        for p in self.parameters:
            val = kwargs[p] if p in kwargs else self.__dict__[p]
            fragment = fragment.replace(p, '(' + str(val) + ')')
        return fragment

    def method(self, **kwargs):
        return ['']

    def survey(self, x0=[0., 0.], alpha0=0.):
        length  = getattr(self, 'length', 0.)
        alpha0 += getattr(self, 'angle',  0.)

        x0[0] += np.cos(alpha0) * length
        x0[1] += np.sin(alpha0) * length

        return x0.copy(), alpha0

class Marker(Element):
    def __init__(self, label):
        self.label = label
        self.parameters = []

    def dump(self, style):
        if style == 'mad' or style == 'elegant':
            return f'{self.label}: MARKER;\n'
        if style == 'at':
            return f"{self.label} = marker('{self.label}', 'IdentityPass');\n"
        if style == 'opa':
            return f'{self.label} : opticsmarker;\r\n'

class Drift(Element):
    def __init__(self, label, length=0.):
        self.label  = label
        self.length = length
        self.parameters = ['length']

    def method(self, flags=0, fracture=[], length=None):
        length = length or self.length

        code = []
        f0 = 0.
        for f in (fracture + [1.]): #Add last point to get the entire pass
            code.append(methods.drift(flags, f'{length} * {(f - f0)}'))
            f0 = f
        return code

    def dump(self, style):
       if style == 'mad' or style == 'elegant':
           return f'{self.label}: DRIFT, L={self.length};\n'
       if style == 'at':
           return f"{self.label} = drift('{self.label}', {self.length}, 'DriftPass');\n"
       if style == 'opa':
           return f"{self.label} : drift, l = {self.length};\r\n"

class Multipole(Element):
    def __init__(self, label, knl=[], ksl=[], dx=0., dy=0.):
        self.label = label
        self.knl   = knl
        self.ksl   = ksl
        self.dx    = dx
        self.dy    = dy
        self.parameters = ['knl', 'ksl', 'dx', 'dy']

    #Multipole of order N can be accessed also as self.kN and self.kNs,
    #this is useful for parametrization (vectors are not implemented)
    def __setattr__(self, name, value):
        if name[0] == 'k' and name[1:].isdigit(): #Normal
            self.knl[int(name[1:])] = value
            return
        if name[0] == 'k' and name[-1] == 's' and name[1:-1].isdigit(): #Skew
            self.ksl[int(name[1:-1])] = value
            return

        super.__setattr__(self, name, value)

    def __getattr__(self, name):
        if name[0] == 'k' and name[1:].isdigit(): #Normal
            return self.knl[int(name[1:])]
        if name[0] == 'k' and name[-1] == 's' and name[1:-1].isdigit(): #Skew
            return self.ksl[int(name[1:-1])]

        return super.__getattr__(self, name)

    def method(self, flags=0, fracture=[], knl=None, ksl=None, dx=None, dy=None, **kwargs):
        if knl == None: knl = self.knl.copy() #Explicit comparison against None
        if ksl == None: ksl = self.ksl.copy() #to discriminate from []
        dx  = dx  or self.dx
        dy  = dx  or self.dy

        for prm, value in kwargs.items():
            if prm[0] == 'k' and prm[1:].isdigit(): #Normal
                knl[int(prm[1:])] = value
            if prm[0] == 'k' and prm[-1] == 's' and prm[1:-1].isdigit(): #Skew
                ksl[int(prm[1:-1])] = value

        code  = methods.align(dx, dy)
        code += methods.kick(flags, knl, ksl)
        code += methods.align(f'-({dx})', f'-({dy})')
        return [code]

    def dump(self, style):
       if style == 'mad' or style == 'elegant':
           knl = ', '.join([str(k) for k in self.knl])
           ksl = ', '.join([str(k) for k in self.ksl])
           return f'{self.label}: MULTIPOLE, KNL={{ {knl} }}, KNL={{ {ksl} }};\n'

       if style == 'at' or style == 'opa':
           raise Exception('Export Multipoles not supported for at and opa')

class Quadrupole(Element):
    def __init__(self, label, slices=None, length=0., k1=0., dx=0., dy=0.):
        self.label  = label
        self.slices = slices if slices else ufo.DEFAULT_QUADRUPOLE_SLICES
        self.length = length
        self.k1     = k1
        self.dx     = dx
        self.dy     = dy
        self.parameters = ['length', 'k1', 'dx', 'dy']

    def method(self, flags=0, fracture=[], slices=None, length=None, k1=None, dx=None, dy=None):
        slices = slices or self.slices
        length = length or self.length
        k1     = k1     or self.k1
        dx     = dx     or self.dx
        dy     = dx     or self.dy
 
        code = []
        f0 = 0.
        for f in (fracture + [1.]): #Add last point to get the entire pass
            code.append(methods.quadrupole(flags, slices, f'{length} * {(f - f0)}', k1))
            f0 = f

        code[0]  = methods.align(dx, dy) + code[0]
        code[-1] = code[-1] + methods.align(f'-({dx})', f'-({dy})')
        return code

    def dump(self, style):
        if style == 'mad' or style == 'elegant':
            return f'{self.label}: QUADRUPOLE, L={self.length}, K1={self.k1};\n'
        if style == 'at':
            return f"{self.label} = quadrupole('{self.label}', {self.length}, {self.k1},'StrMPoleSymplectic4Pass');\n"
        if style == 'opa':
            return f"{self.label} : quadrupole, l = {self.length}, k = {self.k1};\r\n"

class Sbend(Element):
    def __init__(self, label, slices=None, length=0., angle=0., k1=0., e1=0., e2=0., dx=0., dy=0.):
        self.label  = label
        self.slices = slices if slices else ufo.DEFAULT_BEND_SLICES
        self.length = length
        self.angle  = angle
        self.k1     = k1
        self.e1     = e1
        self.e2     = e2
        self.dx     = dx
        self.dy     = dy
        self.parameters = ['length', 'angle', 'k1', 'e1', 'e2', 'dx', 'dy']

    def method(self, flags=0, fracture=[], slices=None, length=None, angle=None, k1=None, e1=None, e2=None, dx=None, dy=None):
        length = length or self.length
        slices = slices or self.slices
        angle  = angle  or self.angle
        k1     = k1     or self.k1
        e1     = e1     or self.e1
        e2     = e2     or self.e2
        dx     = dx     or self.dx
        dy     = dx     or self.dy
 
        code = []
        f0 = 0.
        for f in (fracture + [1.]): #Add last point to get the entire pass
            code.append(methods.sbend(flags, slices, f'{length} * {(f - f0)}', f'{angle} * {(f - f0)}', k1))
            f0 = f

        code[0]   = methods.align(dx, dy) + methods.edge(flags, length, angle, e1) + code[0]
        code[-1] += methods.edge(flags, length, angle, e2) + methods.align(f'-({dx})', f'-({dy})')
        return code

    def dump(self, style):
        if style == 'mad' or style == 'elegant':
            return f'{self.label}: SBEND, L={self.length}, ANGLE={self.angle}, K1={self.k1}, E1={self.e1}, E2={self.e2};\n'
        if style == 'at':
            print('NOT SUPPORTED')
        if style == 'opa':
            print('NOT SUPPORTED')

class Rbend(Element):
    def __init__(self, label, slices=None, length=0., angle=0., k1=0., e1=0., e2=0., dx=0., dy=0.):
        self.label  = label
        self.slices = slices if slices else ufo.DEFAULT_BEND_SLICES
        self.length = length
        self.angle  = angle
        self.k1     = k1
        self.e1     = e1
        self.e2     = e2
        self.dx     = dx
        self.dy     = dy
        self.parameters = ['length', 'angle', 'k1', 'e1', 'e2', 'dx', 'dy']

    def method(self, flags=0, fracture=[], slices=None, length=None, angle=None, k1=None, e1=None, e2=None, dx=None, dy=None):
        length = length or self.length
        slices = slices or self.slices
        angle  = angle  or self.angle
        k1     = k1     or self.k1
        e1     = e1     or self.e1
        e2     = e2     or self.e2
        dx     = dx     or self.dx
        dy     = dx     or self.dy
 
        code = []
        f0 = 0.
        for f in (fracture + [1.]): #Add last point to get the entire pass
            code.append(methods.sbend(flags, slices, f'{length} * {(f - f0)}', f'{angle} * {(f - f0)}', k1))
            f0 = f

        code[0]   = methods.align(dx, dy) + methods.edge(flags, length, angle, f'{e1} + ({angle}*0.5)') + code[0]
        code[-1] += methods.edge(flags, length, angle, f'{e2} + ({angle}*0.5)') + methods.align(f'-({dx})', f'-({dy})')

        return code

    def dump(self, style):
        if style == 'mad' or style == 'elegant':
            return f'{self.label}: RBEND, L={self.length}, ANGLE={self.angle}, K1={self.k1};\n'
        if style == 'at':
            return f"{self.label} = rbend('{self.label}', {self.length}, {self.angle}, 0, 0, {self.k1}, 'BndMPoleSymplectic4Pass');\n"
            return f"{self.label} = rbend('{self.label}', {self.length}, {self.angle}, 0, 0, {self.k1}, 'BndMPoleSymplectic4Pass');\n"
        if style == 'opa':
            print('ERROR RBEND NOT iMPLEMENTED IN OPA')
            exit()
            t  = self.ANGLE / 2. / np.pi * 360.
            t1 = self.e1 / 2. / np.pi * 360.
            t2 = self.e2 / 2. / np.pi * 360.
            return f"{self.label} : bending, l = {self.length}, t = {t}, k = {self.k1}, t1 = {t1}, t2 = {t2};\r\n"

class Sextupole(Element):
    def __init__(self, label, length=0., slices=None, k2=0., k2s=0., dx=0., dy=0.):
        self.label  = label
        self.length = length
        self.slices = slices if slices else ufo.DEFAULT_SEXTUPOLE_SLICES
        self.k2     = k2
        self.k2s    = k2s
        self.dx     = dx
        self.dy     = dy
        self.parameters = ['length', 'k2', 'k2s', 'dx', 'dy']

    def method(self, flags=0, fracture=[], length=None, slices=None, k2=None, k2s=None, dx=None, dy=None):
        length = length or self.length
        slices = slices or self.slices
        k2     = k2     or self.k2
        k2s    = k2s    or self.k2s
        dx     = dx     or self.dx
        dy     = dx     or self.dy

        code = []
        f0 = 0.
        for f in (fracture + [1.]): #Add last point to get the entire pass
            code.append(methods.teapot(flags, f'{length} * {(f - f0)}', slices, [0., 0., k2], [0., 0., k2s]))
            f0 = f

        code[0]  = methods.align(dx, dy) + code[0]
        code[-1] = code[-1] + methods.align(f'-({dx})', f'-({dy})')
        return code

    def dump(self, style):
        if style == 'mad' or style == 'elegant':
            return f'{self.label}: SEXTUPOLE, L={self.length}, K2={self.k2};\n'
        if style == 'at':
            return f"{self.label} = sextupole('{label}', {self.length}, {self.k2},'StrMPoleSymplectic4Pass');\n"
        if style == 'opa':
            return f"{self.label} : sextupole, l = {self.length}, k = {self.k2 * .5};\r\n"

class Octupole(Element):
    def __init__(self, label, length=0., slices=None, k3=0., k3s=0., dx=0., dy=0.):
        self.label  = label
        self.length = length
        self.slices = slices if slices else ufo.DEFAULT_SEXTUPOLE_SLICES
        self.k3     = k3
        self.k3s    = k3s
        self.dx     = dx
        self.dy     = dy
        self.parameters = ['length', 'k3', 'k3s', 'dx', 'dy']

    def method(self, flags=0, fracture=[], length=None, slices=None, k3=None, k3s=None, dx=None, dy=None):
        length = length or self.length
        slices = slices or self.slices
        k3     = k3     or self.k3
        k3s    = k3s    or self.k3s
        dx     = dx     or self.dx
        dy     = dx     or self.dy

        code = []
        f0 = 0.
        for f in (fracture + [1.]): #Add last point to get the entire pass
            code.append(methods.teapot(flags, f'{length} * {(f - f0)}', slices, [0., 0., 0., k3], [0., 0., 0., k3s]))
            f0 = f

        code[0]  = methods.align(dx, dy) + code[0]
        code[-1] = code[-1] + methods.align(f'-({dx})', f'-({dy})')
        return code

    def dump(self, style):
        if style == 'mad' or style == 'elegant':
            return f'{self.label}: OCTUPOLE, L={self.length}, K3={self.k3};\n'
        if style == 'at':
            return f"{self.label} = octupole('{label}', {self.length}, {self.k3},'StrMPoleSymplectic4Pass');\n"
        if style == 'opa':
            return f"{self.label} : octupole, l = {self.length}, k = {self.k3 * .5};\r\n"

class Lattice():
    def __init__(self, path=None):
        if path:
            self.__dict__.update(mad.load(path))

    def __repr__(self):
        return '\n'.join([f'{repr(self.__dict__[p])}' for p in self.__dict__])

class Line(list):
    def __init__(self, label):
        self.label = label #Same label convention as elements
        super().__init__()

    def __str__(self):
        return self.label 

    def __repr__(self):
        return f'{self.label: <8}: LINE         {[str(element) for element in self]}'

    #Thess are implemented as properties to match elements behaviour
    def __get_count(self):
        return sum([getattr(e, 'count', 1.) for e in self])

    def __get_length(self):
        return sum([getattr(e, 'length', 0.) for e in self])

    def __get_angle(self):
        return sum([getattr(e, 'angle', 0.) for e in self])

    count  = property(__get_count,  None)
    length = property(__get_length, None)
    angle  = property(__get_angle,  None)

    def find(self, what):
        where=[]
        for idx, element in enumerate(self.flatten()):
            if what(element):
                where.append(idx)

        return where

    def locate(self, element): #Warning maybe this could be implemented in a smarter way?
        if element == -1: return self.length
        flat_line = self.flatten()
        s = 0
        for idx, e in enumerate(flat_line):
            if int(element) == idx:
                s += (element - int(element)) * e.length
                return s
            s += getattr(e, 'length', 0.)

    def survey(self, x0=[0., 0.], alpha0=0.):
        x = []
        alpha = []
        for e in self:
            x0, alpha0 = e.survey(x0=x0, alpha0=alpha0)
            x.append(x0)
            alpha.append(alpha0)

        return x, alpha    
        #return [x = e.survey(x=x, alpha=alpha) for e in self]

    def method(self, fracture=[], flags=0, parameters=[]):
        line = self.flatten()
        code = ['']
        for idx, element in enumerate(line):
            #Check where to fracture
            where = [f - idx for f in fracture if f != '*' and idx <= f < idx+1]
            if 0 in where: #Produces a fracture point at begin of element
                where.remove(0)
                code.append('') #Dummy code to induce a fracture point

            #Current element parameters, 'i' is the absolute parameter index
            #Match elemet label or index
            prms = [i for (i, p) in enumerate(parameters) if p[0]==element.label or p[0]==idx]
            prms = {parameters[i][1]: f'parameters[idx][{i}]' for i in prms}

            prms['flags']    = flags #also flags and fracture need to 
            prms['fracture'] = where #be passed to the element method()

            fragment = element.method(**prms)
            code[-1] += fragment.pop(0)
            code += fragment

        return code

    def flatten(self):
        flat = [e.flatten() if getattr(e, 'flatten', None) else [e] for e in self]
        return sum(flat, []) #Flatten list of lists

    def dump(self, style):
        if style == 'mad' or style == 'elegant':
            line = f'{self.label}: LINE=('
        if style == 'at':
            line = f'{self.label} = ['
        if style == 'opa':
            line = f'{self.label} : '
        not_first = 0
        for element in self:
            if not_first:
                    line += ', '
            not_first = 1
            line += f'{element.label}'

        if style == 'mad':
            line += ');\n\n'
        if style == 'elegant':
            line += ')\n\n'
        if style == 'at':
            line += '];\n\n'
        if style == 'opa':
            line += ';\r\n'

        return line

class Beam():
    def __init__(self, energy=3e9, particle_mass=electron_mass, bunch_charge=1e-9, beam_current=0.25,
                 ex=1e-9, ey=1e-9, bunch_length=6e-3, energy_spread=1e-3):
        self.energy = energy
        self.particle_mass   = particle_mass
        self.gamma  = (1. + energy) / particle_mass
        self.bunch_charge = bunch_charge
        self.beam_current = beam_current
        self.ex = ex
        self.ey = ey
        self.bunch_length = bunch_length
        self.energy_spread = energy_spread

def dump(lattice, path, style='mad', beam=Beam()): #style can be 'mad', 'elegant', 'at' or opa
    dump_file = open(path, 'w')
    if style == 'at':
        dump_file.write('function {}\n\n'.format(path.split(".")[0]))
        dump_file.write('global FAMLIST THERING GLOBVAL;\n')
        dump_file.write(f'E_0 = {beam.mass};\n')
        dump_file.write(f'GLOBVAL.E0 = {beam.energy}-E_0;\n')
        dump_file.write("GLOBVAL.LatticeFile = '{}';\n".format(path.split(".")[0]))
        dump_file.write('FAMLIST = cell(0);\n')
        dump_file.write('THERING = cell(0);\n\n')

    if style == 'opa':
        dump_file.write(f'energy = {beam.energy * 1e-9};\r\n')

    for element in lattice.__dict__.values():
        dump_file.write(element.dump(style))
    
    if style == 'opa':
        dump_file.write('\r')
    dump_file.write('\n')

    if style == 'at':
        dump_file.write('\n\nbuildlat(RING);\n')
        dump_file.write('% Set all magnets to same energy\n')
        dump_file.write("THERING = setcellstruct(THERING, 'Energy', 1:length(THERING), GLOBVAL.E0);\n")
        dump_file.write("evalin('caller', 'global THERING FAMLIST GLOBVAL');\n")
        dump_file.write('if nargout\n    varargout{1} = THERING;\nend\nend')

    dump_file.close()

