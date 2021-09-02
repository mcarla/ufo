#   You can redistribute this file and/or modify it under the terms of the GNU
#   General Public License GPLv3 (or later), as published by the Free Software
#   Foundation. This file is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.

#   Authors: M. Carla'

LINEAR           = 1 << 0
FIVED            = 1 << 1
EXACT            = 1 << 2
KICK             = 1 << 3
RADIATION        = 1 << 4
DOUBLE_PRECISION = 1 << 5
ACHROMATIC       = 1 << 6

def align(dx, dy): #MAD-X definitions
    fragment  = (f'x -= {dx};\n'
                 f'y -= {dy};\n')
 
    return fragment

def drift(flags, length):
    if flags & EXACT:
        fragment  = ('{\n'
                     '    ufloat opdp, l;\n'

                     '    opdp = 1. + dp;\n'
                    f'    l = {length} / sqrt(opdp*opdp - px*px - py*py);\n'
                     '    x += px * l;\n'
                     '    y += py * l;\n'
                     '}\n')
 
    else:
        fragment = (f'x += px * {(length)};\n'
                    f'y += py * {(length)};\n')

    return fragment

def kick(flags, knl, ksl):
    fragment = ''

    if knl:
        n = len(knl) #Avoid non-linear terms when LINEAR
        max_order = min(n, 2) if (flags & LINEAR) else n 

        dpx = f'{knl[max_order - 1]} {"" if flags & ACHROMATIC else "* oodppo"}'
        dpy =  '0.'

        for order in range(max_order - 1, 0, -1):
            aux = f'({dpx} * x - ({dpy} * y)) / {order}'
            dpy = f'({dpx} * y + ({dpy} * x)) / {order}'
            dpx = f'({knl[order - 1]} {"" if flags & ACHROMATIC else "* oodppo"} + {aux})'

        fragment = f'px -= {dpx};\npy += {dpy};\n\n'

    if ksl:
        n = len(ksl) #Avoid non-linear terms when LINEAR
        max_order = min(n, 2) if (flags & LINEAR) else n 

        dpx =  '0.'
        dpy = f'{ksl[max_order - 1]} {"" if flags & ACHROMATIC else "* oodppo"}'

        for order in range(max_order - 1, 0, -1):
            aux = f'({dpx} * y + ({dpy} * x)) / {order}'
            dpx = f'({dpx} * x - ({dpy} * y)) / {order}'
            dpy = f'({ksl[order - 1]} {"" if flags & ACHROMATIC else "* oodppo"} + {aux})'
        fragment += f'px -= {dpx};\npy += {dpy};\n\n'

    return fragment

#https://accelconf.web.cern.ch/IPAC2013/papers/mopwo027.pdf
def teapot(flags, length, slices, knl, ksl): #Using Teapot slicing
    edge  = f'({length} / 2.)' if slices == 1 else f'({length} / (2 * (1 + {slices})))' #Edge drifts length
    inner = f'({length} - 2. * {edge}) / {slices - 1.}' #Inner drifts length

    weak_edge  = f'px -= x * (({knl[0]}) * ({knl[0]}) * ({edge})  / ({length}) / ({length}));\n' #Weak focusing
    weak_inner = f'px -= x * (({knl[0]}) * ({knl[0]}) * ({inner}) / ({length}) / ({length}));\n'

    knl = [f'({k}*{length}/{slices})' for k in knl]
    ksl = [f'({k}*{length}/{slices})' for k in ksl]

    fragment  = weak_edge + drift(flags, edge)
    for count in range(slices - 1):
        fragment += kick(flags, knl, ksl)
        fragment += weak_inner
        fragment += drift(flags, inner)
    fragment += kick(flags, knl, ksl) + drift(flags, edge) + weak_edge

    return fragment

def quadrupole(flags, slices, length, k1):
    if flags & KICK:
        fragment = teapot(flags, length, slices, [0., k1], [])

    else:
        fragment = ('{\n'
                    '    ufloat x0  = x;\n'
                    '    ufloat px0 = px;\n'
                    '    ufloat y0  = y;\n'
                    '    ufloat py0 = py;\n'

                   f'    ufloat k         = {k1};\n'
                   f'    ufloat k2        = sqrt(fabs(k {"" if flags & ACHROMATIC else "* oodppo"}));\n'
                   f'    ufloat k2l       = {length} * k2;\n'
                    '    ufloat C         = cos( k2l);\n'
                    '    ufloat CH        = cosh(k2l);\n'
                    '    ufloat S         = sin( k2l);\n'
                    '    ufloat SH        = sinh(k2l);\n'

                    '    if (k > 0.) {\n'
                    '        x  =  C  * x0      + S  * px0 / k2;\n'
                    '        px = -S  * x0 * k2 + C  * px0;\n'
                    '        y  =  CH * y0      + SH * py0 / k2;\n'
                    '        py =  SH * y0 * k2 + CH * py0;\n'
                    '    } else {\n'
                    '        x  =  CH * x0      + SH * px0 / k2;\n'
                    '        px =  SH * x0 * k2 + CH * px0;\n'
                    '        y  =  C  * y0      + S  * py0 / k2;\n'
                    '        py = -S  * y0 * k2 + C  * py0;\n'
                    '   }\n'
                    '}\n')

    return fragment

def sbend(flags, slices, length, angle, k1):
    if flags & KICK:
        fragment = teapot(flags, length, slices, [angle, k1], [])
    else:
        fragment = ('{\n'
                    '    ufloat x0  = x;\n'
                    '    ufloat px0 = px;\n'
                    '    ufloat y0  = y;\n'
                    '    ufloat py0 = py;\n'

                   f'    ufloat curvature = {angle} / ({length}) {"" if flags & ACHROMATIC else "* oodppo"};\n'
                   f'    ufloat k         = {k1} {"" if flags & ACHROMATIC else "* oodppo"};\n'
                    '    ufloat k2        = sqrt(fabs(k));\n'
                   f'    ufloat k2l       = {length} * k2;\n'
                    '    ufloat C         = cos( k2l);\n'
                    '    ufloat CH        = cosh(k2l);\n'
                    '    ufloat S         = sin( k2l);\n'
                    '    ufloat SH        = sinh(k2l);\n'

                    '    if (k > 0.) {\n'
                    '        y  =  CH * y0      + SH * py0 / k2;\n'
                    '        py =  SH * y0 * k2 + CH * py0;\n'
                    '    }\n'

                    '    if (k < 0.) {\n'
                    '        y  =  C  * y0      + S  * py0 / k2;\n'
                    '        py = -S  * y0 * k2 + C  * py0;\n'
                    '    }\n'

                    '    if (k == 0.)\n'
                   f'        y += {length} * py0;\n'

                   f'    k += curvature * curvature {"" if flags & ACHROMATIC else "/ oodppo"};\n' #Horizontal plane includes weak focusing 
                    '    k2  = sqrt(fabs(k));\n'
                   f'    k2l = {length} * k2;\n'
                    '    C   = cos( k2l);\n'
                    '    CH  = cosh(k2l);\n'
                    '    S   = sin( k2l);\n'
                    '    SH  = sinh(k2l);\n'

                    '    if (k > 0.) {\n'
                    '        x   =  C  * x0      + S  * px0 / k2;\n'
                   f'        x  += dp * curvature * (1. - C) / fabs(k);\n'
                    '        px  = -S  * x0 * k2 + C  * px0;\n'
                   f'        px += dp * curvature * S / k2;\n'
                    '    }\n'

                    '    if (k < 0.) {\n'
                    '        x  =  CH * x0      + SH * px0 / k2;\n'
                   f'        x +=  dp * curvature * (CH - 1.) / fabs(k);\n'
                    '        px =  SH * x0 * k2 + CH * px0;\n'
                   f'        px += dp * curvature * SH / k2;\n'
                    '    }\n'

                    '    if (k == 0.)\n'
                   f'        x += {length} * px0;\n'
                    '}\n')

    return fragment

def edge(flags, length, angle, e):
    fragment = ('{\n'
               f'    ufloat curvature = {angle} / ({length}) {"" if flags & ACHROMATIC else "* oodppo"};\n'

               f'    px += x * curvature * tan({e});\n'
               f'    py -= y * curvature * tan({e});\n'
                '}\n')

    return fragment

#def cavity(flags, field, frequency, lag):
#    #field is normalized by pc (reference beam energy)
#    #frequency is in Hz
#    omega = '({frequency}) / {speed_of_light}'
#    return f'{dp += field * sin({2. * np.pi} * ({lag} - ({omega}) * z))};\n'
