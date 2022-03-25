#   You can redistribute this file and/or modify it under the terms of the GNU
#   General Public License GPLv3 (or later), as published by the Free Software
#   Foundation. This file is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.

#   Authors: M. Carla'

import numpy as np
from . import methods
from . import mad

import ufo

def __dk_from_kwargs__(kwargs, dkn, dks):
    for prm, value in kwargs.items():
        if prm[:2] == 'dk' and prm[2:].isdigit(): #Normal
            idx = int(prm[2:])
            dkn += [0.] * (idx + 1 - len(dkn)) #pad list to length
            dkn[idx] = value
        if prm[:2] == 'dk' and prm[-1] == 's' and prm[2:-1].isdigit(): #Skew
            idx = int(prm[2:-1])
            dks += [0.] * (idx + 1 - len(dks)) #pad list to length
            dks[idx] = value

class Element():
    def __init__(self):
        pass

    def __str__(self):
        return self.label 

    def __repr__(self):  
        msg = f'{self.label: <8}: {type(self).__name__.upper(): <12} '
        return msg + '   '.join([f'{p}={self.__dict__[p]}' for p in self.parameters])

    #Multipolar field error of order N can be accessed also as self.dkN and self.dkNs,
    #this is useful for parametrization (vectors are not implemented)
    def __setattr__(self, name, value):
        if name[:2] == 'dk' and name[2:].isdigit(): #Normal
            self.dkn[int(name[2:])] = value
            return
        if name[:2] == 'dk' and name[-1] == 's' and name[2:-1].isdigit(): #Skew
            self.dks[int(name[2:-1])] = value
            return

        super.__setattr__(self, name, value)

    def __getattr__(self, name):
        if name[:2] == 'dk' and name[2:].isdigit(): #Normal
            return self.dkn[int(name[2:])]
        if name[:2] == 'dk' and name[-1] == 's' and name[2:-1].isdigit(): #Skew
            return self.dks[int(name[2:-1])]

        return super.__getattr__(self, name)

    def method(self, **kwargs):
        return ['']

class Marker(Element):
    """
    A marker element

    label : str -> Name of the element
    """
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
    """
    A drift element. The flag EXACT enables exact hamiltonian pass method

    label  : str   -> Name of the element
    length : float -> Length of the drift
    """
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
    """
    A thin multipole element.

    label    : str   -> Name of the element
    knl, ksl : list  -> List of normal and skew multipoles
    dx, dy   : float -> Horizontal and vertical alignment offset
    """
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
        dx = dx or self.dx
        dy = dy or self.dy

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
    """
    A thick normal quadrupole element.

    label    : str   -> Name of the element
    slices   : float -> Number of slices used for 'tea pot' expansion (if KICK flag is set)
    length   : float -> Length of the element
    k1       : float -> Strength of the quadrupolar component
    dx, dy   : float -> Horizontal and vertical alignment offset
    dkn, dks : list  -> List of normal and skew multipolar field errors
    """
    def __init__(self, label, slices=None, length=0., k1=0., dx=0., dy=0., dkn=[], dks=[]):
        self.label  = label
        self.slices = slices if slices else ufo.DEFAULT_QUADRUPOLE_SLICES
        self.length = length
        self.k1     = k1
        self.dx     = dx
        self.dy     = dy
        self.dkn    = dkn
        self.dks    = dks
        self.parameters = ['length', 'k1', 'dx', 'dy', 'dkn', 'dks']

    def method(self, flags=0, fracture=[], slices=None, length=None, k1=None, dx=None, dy=None, dkn=None, dks=None, **kwargs):
        slices = slices or self.slices
        length = length or self.length
        k1     = k1     or self.k1
        dx     = dx     or self.dx
        dy     = dy     or self.dy
        if dkn == None: dkn = self.dkn.copy() #Explicit comparison against None
        if dks == None: dks = self.dks.copy() #to discriminate from []

        __dk_from_kwargs__(kwargs, dkn, dks)

        code = []
        f0 = 0.
        for f in (fracture + [1.]): #Add last point to get the entire pass
            code.append(methods.quadrupole(flags, slices, f'{length} * {(f - f0)}', k1, dkn, dks))
            f0 = f

        code[ 0] = methods.align(dx, dy) + code[0]
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
    """
    An S type bending magnet with edge focusing.

    label    : str   -> Name of the element
    slices   : float -> Number of slices used for 'tea pot' expansion (if KICK flag is set)
    length   : float -> Length of the element
    angle    : float -> Bending angle in radiants
    e1, e2   : float -> Entry and exit angle
    hgap     : float -> Half gap of the magnet
    fint     : float -> Entrance and exit fringe field integral
    k1       : float -> Strength of the quadrupolar component
    dx, dy   : float -> Horizontal and vertical alignment offset
    dkn, dks : list  -> List of normal and skew multipolar field errors
    """
    def __init__(self, label, slices=None, length=0., angle=0., k1=0., e1=0., e2=0., hgap=0., fint=0., dx=0., dy=0., dkn=[], dks=[]):
        self.label  = label
        self.slices = slices if slices else ufo.DEFAULT_BEND_SLICES
        self.length = length
        self.angle  = angle
        self.k1     = k1
        self.e1     = e1
        self.e2     = e2
        self.hgap   = hgap
        self.fint   = fint
        self.dx     = dx
        self.dy     = dy
        self.dkn    = dkn
        self.dks    = dks
        self.parameters = ['length', 'angle', 'k1', 'e1', 'e2', 'hgap', 'fint', 'dx', 'dy', 'dkn', 'dks']

    def method(self, flags=0, fracture=[], slices=None, length=None, angle=None, k1=None, e1=None, e2=None, hgap=None, fint=None, dx=None, dy=None, dkn=None, dks=None, **kwargs):
        length = length or self.length
        slices = slices or self.slices
        angle  = angle  or self.angle
        k1     = k1     or self.k1
        e1     = e1     or self.e1
        e2     = e2     or self.e2
        hgap   = hgap   or self.hgap
        fint   = fint   or self.fint
        dx     = dx     or self.dx
        dy     = dy     or self.dy
        if dkn == None: dkn = self.dkn.copy() #Explicit comparison against None
        if dks == None: dks = self.dks.copy() #to discriminate from []

        __dk_from_kwargs__(kwargs, dkn, dks)

        code = []
        f0 = 0.
        for f in (fracture + [1.]): #Add last point to get the entire pass
            code.append(methods.sbend(flags, slices, f'{length} * {(f - f0)}', f'{angle} * {(f - f0)}', k1, dkn, dks))
            f0 = f

        fringe = f'({hgap} * {fint})'
        code[ 0]  = methods.align(dx, dy) + methods.edge(flags, length, angle, e1, fringe) + code[0]
        code[-1] += methods.edge(flags, length, angle, e2, fringe) + methods.align(f'-({dx})', f'-({dy})')
        return code

    def dump(self, style):
        if style == 'mad' or style == 'elegant':
            return f'{self.label}: SBEND, L={self.length}, ANGLE={self.angle}, K1={self.k1}, E1={self.e1}, E2={self.e2};\n'
        if style == 'at':
            print('NOT SUPPORTED')
        if style == 'opa':
            print('NOT SUPPORTED')

class Rbend(Element):
    """
    A rectangular bending magnet with edge focusing.
    
    label    : str   -> Name of the element
    slices   : float -> Number of slices used for 'tea pot' expansion (if KICK flag is set)
    length   : float -> Length of the element
    angle    : float -> Bending angle in radiants
    e1, e2   : float -> Entry and exit angle
    hgap     : float -> Half gap of the magnet
    fint     : float -> Entrance and exit fringe field integral
    k1       : float -> Strength of the quadrupolar component
    dx, dy   : float -> Horizontal and vertical alignment offset
    dkn, dks : list  -> List of normal and skew multipolar field errors
    """
    def __init__(self, label, slices=None, length=0., angle=0., k1=0., e1=0., e2=0., hgap=0, fint=0, dx=0., dy=0., dkn=[], dks=[]):
        self.label  = label
        self.slices = slices if slices else ufo.DEFAULT_BEND_SLICES
        self.length = length
        self.angle  = angle
        self.k1     = k1
        self.e1     = e1
        self.e2     = e2
        self.hgap   = hgap
        self.fint   = fint
        self.dx     = dx
        self.dy     = dy
        self.dkn    = dkn
        self.dks    = dks
        self.parameters = ['length', 'angle', 'k1', 'e1', 'e2', 'hgap', 'fint', 'dx', 'dy', 'dkn', 'dks']

    def method(self, flags=0, fracture=[], slices=None, length=None, angle=None, k1=None, e1=None, e2=None, hgap=None, fint=None, dx=None, dy=None, dkn=None, dks=None, **kwargs):
        length = length or self.length
        slices = slices or self.slices
        angle  = angle  or self.angle
        k1     = k1     or self.k1
        e1     = e1     or self.e1
        e2     = e2     or self.e2
        hgap   = hgap   or self.hgap
        fint   = fint   or self.fint
        dx     = dx     or self.dx
        dy     = dy     or self.dy
        if dkn == None: dkn = self.dkn.copy() #Explicit comparison against None
        if dks == None: dks = self.dks.copy() #to discriminate from []

        __dk_from_kwargs__(kwargs, dkn, dks)
 
        code = []
        f0 = 0.
        for f in (fracture + [1.]): #Add last point to get the entire pass
            code.append(methods.sbend(flags, slices, f'{length} * {(f - f0)}', f'{angle} * {(f - f0)}', k1, dkn, dks))
            f0 = f

        fringe = f'({hgap} * {fint})'
        code[0]   = methods.align(dx, dy) + methods.edge(flags, length, angle, f'{e1} + ({angle}*0.5)', fringe) + code[0]
        code[-1] += methods.edge(flags, length, angle, f'{e2} + ({angle}*0.5)', fringe) + methods.align(f'-({dx})', f'-({dy})')

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
    """
    A thick sextupole ('Tea-pot expansion used)
 
    label    : str   ->  Name of the element
    slices   : float ->  Number of slices used for 'tea pot' expansion
    length   : float ->  Length of the element
    k2       : float ->  Normal sextupolar component
    k2s      : float ->  Skew sextupolar component
    dx, dy   : float ->  Horizontal and vertical alignment offset
    dkn, dks : list  -> List of normal and skew multipolar field errors
    """
    def __init__(self, label, length=0., slices=None, k2=0., k2s=0., dx=0., dy=0., dkn=[], dks=[]):
        self.label  = label
        self.length = length
        self.slices = slices if slices else ufo.DEFAULT_SEXTUPOLE_SLICES
        self.k2     = k2
        self.k2s    = k2s
        self.dx     = dx
        self.dy     = dy
        self.dkn    = dkn
        self.dks    = dks
        self.parameters = ['length', 'k2', 'k2s', 'dx', 'dy', 'dkn', 'dks']

    def method(self, flags=0, fracture=[], length=None, slices=None, k2=None, k2s=None, dx=None, dy=None, dkn=None, dks=None, **kwargs):
        length = length or self.length
        slices = slices or self.slices
        k2     = k2     or self.k2
        k2s    = k2s    or self.k2s
        dx     = dx     or self.dx
        dy     = dy     or self.dy
        if dkn == None: dkn = self.dkn.copy() #Explicit comparison against None
        if dks == None: dks = self.dks.copy() #to discriminate from []

        __dk_from_kwargs__(kwargs, dkn, dks)

        code = []
        f0 = 0.
        for f in (fracture + [1.]): #Add last point to get the entire pass
            #code.append(methods.teapot(flags, f'{length} * {(f - f0)}', slices, [0., 0., k2], [0., 0., k2s]))
            code.append(methods.sextupole(flags, slices, f'{length} * {(f - f0)}', k2, k2s, dkn, dks))
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
    """
    A thick octupole ('Tea-pot expansion used)
 
    label    : str   -> Name of the element
    slices   : float -> Number of slices used for 'tea pot' expansion
    length   : float -> Length of the element
    k2       : float -> Normal octupolar component
    k2s      : float -> Skew octupolar component
    dx, dy   : float -> Horizontal and vertical alignment offset
    dkn, dks : list  -> List of normal and skew multipolar field errors
    """
    def __init__(self, label, length=0., slices=None, k3=0., k3s=0., dx=0., dy=0., dkn=[], dks=[]):
        self.label  = label
        self.length = length
        self.slices = slices if slices else ufo.DEFAULT_OCTUPOLE_SLICES
        self.k3     = k3
        self.k3s    = k3s
        self.dx     = dx
        self.dy     = dy
        self.dkn    = dkn
        self.dks    = dks
        self.parameters = ['length', 'k3', 'k3s', 'dx', 'dy', 'dkn', 'dks']

    def method(self, flags=0, fracture=[], length=None, slices=None, k3=None, k3s=None, dx=None, dy=None, dkn=None, dks=None, **kwargs):
        length = length or self.length
        slices = slices or self.slices
        k3     = k3     or self.k3
        k3s    = k3s    or self.k3s
        dx     = dx     or self.dx
        dy     = dy     or self.dy
        if dkn == None: dkn = self.dkn.copy() #Explicit comparison against None
        if dks == None: dks = self.dks.copy() #to discriminate from []

        __dk_from_kwargs__(kwargs, dkn, dks)

        code = []
        f0 = 0.
        for f in (fracture + [1.]): #Add last point to get the entire pass
            #code.append(methods.teapot(flags, f'{length} * {(f - f0)}', slices, [0., 0., 0., k3], [0., 0., 0., k3s]))
            code.append(methods.octupole(flags, slices, f'{length} * {(f - f0)}', k3, k3s, dkn, dks))
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

class Wire(Element):
    """
    Kick produced by a wire parallel to the beam.
    The magnitude of the produced kick is k / r.

    label : str   -> Name of the element
    x, y  : float -> Horizontal and vertical position of the wire 
    k     : float -> Field strength mu0 * I / (2pi), with I the wire current
    """
    def __init__(self, label, x=0., y=0., k=0.):
        self.label = label
        self.x     = x
        self.y     = y
        self.k     = k
        self.parameters = ['x', 'y', 'k']

    def method(self, flags=0, fracture=[], x=None, y=None, k=None, **kwargs):
        x    = x     or self.x
        y    = y     or self.y
        k    = k     or self.k

        code = methods.wire(x, y, k)
        return [code] #[methods.wire(x, y, k)]

    def dump(self, style):
        print('Warning: dumping of Wire elements is not supported!')

