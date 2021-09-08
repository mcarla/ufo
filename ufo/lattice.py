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

class Lattice():
    """A class containing elements and lines"""
    def __init__(self, path=None):
        if path:
            self.__dict__.update(mad.load(path))

    def __repr__(self):
        return '\n'.join([f'{repr(self.__dict__[p])}' for p in self.__dict__])

class Line(list):
    """A class derivated from the list class, containing a list of elements or other lines"""
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
        """Return a list elements positons such that 'what(element) == True'"""
        where=[]
        for idx, element in enumerate(self.flatten()):
            if what(element):
                where.append(idx)

        return where

    def locate(self, element): #Warning maybe this could be implemented in a smarter way?
        """Return the position in meters of 'element'"""
        if element == -1: return self.length
        flat_line = self.flatten()
        s = 0
        for idx, e in enumerate(flat_line):
            if int(element) == idx:
                s += (element - int(element)) * e.length
                return s
            s += getattr(e, 'length', 0.)

    def survey(self, x0=[0., 0.], alpha0=0.):
        """Not fully implemented"""
        x = []
        alpha = []
        for e in self:
            x0, alpha0 = e.survey(x0=x0, alpha0=alpha0)
            x.append(x0)
            alpha.append(alpha0)

        return x, alpha    
        #return [x = e.survey(x=x, alpha=alpha) for e in self]

    def method(self, fracture=[], flags=0, parameters=[]):
        """
        Return the OpenCL pass method

        Parameters
        ----------
            fracture : list
                When fracture is provided the returned method will be splitted
                in multiples methods, each relative to a fraction of the line
            flags : int
                Combination of: LINEAR, FIVED, EXACT, KICK, RADIATION,
                                DOUBLE_PRECISION, ACHROMATIC
            parameters : list
                A list of parameters that will not be hardcoded in the method
        """
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
        """Return a flattened version of the line (all lines are expanded into its composing elements)"""
        flat = [e.flatten() if getattr(e, 'flatten', None) else [e] for e in self]
        return sum(flat, []) #Flatten list of lists

    def dump(self, style):
        """Return a text representation of the line compatible with one of mad, elegant, opa or at file formats"""
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
    """Class containg some beam informations"""
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
    """
    Save an antire lattice to a file compatible with mad, elegant, accelerator toolbox or opa

    Parameters
    ----------
        lattice : Lattice
            The lattice to be saved
        path : str
            Path to the saved file
        style : str
            One among: 'mad', 'elegant', 'at', 'opa'
    """
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

