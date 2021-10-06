#   You can redistribute this file and/or modify it under the terms of the GNU
#   General Public License GPLv3 (or later), as published by the Free Software
#   Foundation. This file is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.

#   Authors: M. Carla'

import collections
import re
import ufo

Token = collections.namedtuple('Token', ['type', 'value'])

def tokenize(mad_file):
    token_specification = [
                             #Match scalar   | match vector {a, b, c...}
        ('NUMERIC',         r'(-?\d+(\.\d*)?)|(\s*\{(-?\d+(\.\d*)?\s*,?\s*)*\})'),
        ('OPEN',            r'\('),                     #
        ('CLOSE',           r'\)'),                     #
        ('COLON',           r':'),                      #
        ('COMMA',           r','),                      #
        ('ASSIGN',          r'='),                      # Assignment operator
        ('DEFERRED_ASSIGN', r':='),                     # Assignment operator
        ('END',             r';'),                      # Statement terminator
        ('ID',              r'[A-Za-z_\-][A-Za-z0-9_]*'), # Identifiers
        #('OP',              r'[+\-*/]'),                # Arithmetic operators
        ('NEWLINE',         r'\n'),                     # Line endings
        ('SKIP',            r'[ \t]+'),                 # Skip over spaces and tabs
        ('MISMATCH',        r'.'),                      # Any other character
    ]

    tok_regex = '|'.join('(?P<%s>%s)' % pair for pair in token_specification)
    for mo in re.finditer(tok_regex, mad_file.read()):
        kind = mo.lastgroup
        value = mo.group()

        if kind == 'NUMERIC':
            if '{' in value: #Is a vector
                value = value.replace('{', '').replace('}', '').split(',')
                value = [float(n) for n in value]
            else: #Is a scalar
                value = float(value) #if '.' in value else int(value)
        elif kind == 'NEWLINE':
            continue
        elif kind == 'SKIP':
            continue
        elif kind == 'MISMATCH':
            raise RuntimeError(f'{value!r} unexpected')

        yield Token(kind, value)

def parse_element(tokenizer):
    parameters = {}
    while True:
        token = next(tokenizer)
        if token.type == 'END':
            break
        if token.type != 'COMMA':
            print('There shoud have been a comma, but I forgive you')

        token = next(tokenizer)
        parameter_name = token.value

        token = pop_or_die(tokenizer, 'ASSIGN')
        token = pop_or_die(tokenizer, 'NUMERIC')
        parameters[parameter_name] = token.value

    return parameters

def parse_line(tokenizer, lattice, line):
    token = pop_or_die(tokenizer, 'ASSIGN')
    token = pop_or_die(tokenizer, 'OPEN')

    while True:
            token = pop_or_die(tokenizer, 'ID')
            line.append(lattice[token.value])

            token = next(tokenizer)
            if token.type == 'CLOSE':
                break
            elif token.type != 'COMMA':
                raise RuntimeError(f'{token.value!r} unexpected')

    pop_or_die(tokenizer, 'END')
    return line

def parse_unknown(tokenizer):
    while True:
            token = next(tokenizer)
            if token.type == 'END':
                    break

def pop_or_die(tokenizer, kind):
    token = next(tokenizer)
    if token.type != kind:
            raise RuntimeError(f'{token.value!r} unexpected')
    return token

def load(path, flatten=True):

    elements = {
        'MARKER'    : lambda label:                        ufo.Marker(label),
        'DRIFT'     : lambda label, L=0.:                  ufo.Drift(label, length=L),
        'QUADRUPOLE': lambda label, L=0., K1=0.:           ufo.Quadrupole(label, length=L, k1=K1),
        'SBEND'     : lambda label, L=0., ANGLE=0., K1=0., E1=0., E2=0., HGAP=0., FINT=0.: ufo.Sbend(label, length=L, angle=ANGLE, k1=K1, e1=E1, e2=E2, hgap=HGAP, fint=FINT),
        'RBEND'     : lambda label, L=0., ANGLE=0., K1=0., E1=0., E2=0.: ufo.Rbend(label, length=L, angle=ANGLE, k1=K1, e1=E1, e2=E2),
        'SEXTUPOLE' : lambda label, L=0., K2=0., K2S=0.:   ufo.Sextupole(label, length=L, k2=K2, k2s=K2S),
        'OCTUPOLE'  : lambda label, L=0., K3=0., K3S=0.:   ufo.Octupole(label, length=L, k3=K3, k3s=K3S),
        'MULTIPOLE' : lambda label, KNL=[], KSL=[]:        ufo.Multipole(label, knl=KNL, ksl=KSL)
    }

    try:
        mad_file = open(path, 'r')
    except:
        print('Error: could not open input file!')

    lattice = {}
    tokenizer = tokenize(mad_file)
    while True:
        try:
            token = next(tokenizer)
        except StopIteration:
            break

        if token.type != 'ID':
            raise RuntimeError(f'{token.value!r} unexpected')

        label = token.value

        token = pop_or_die(tokenizer, 'COLON')
        token = pop_or_die(tokenizer, 'ID')
        kind = token.value.upper()
        if kind in elements:
            parameters = parse_element(tokenizer)
            lattice[label] = elements[kind](label, **parameters)
        elif kind == 'LINE':
            line = parse_line(tokenizer, lattice, ufo.lattice.Line(label))
            lattice[label] = line
        else:
            print(f'Warning: unknown statement \'{kind}\'')
            parse_unknown(tokenizer)

    mad_file.close()
    return lattice

