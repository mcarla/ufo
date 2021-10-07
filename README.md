                            U.F.O.
 
 
                          \__/  ^__^
                          (oo)  (oo)
                         _/--\  /--\_
                   _.--===0=0====0=0===--._
                  (________________________)
                       /  \________/  \
                      /                \
 
 
       An Unreliable, but (Undoubtedly) Fast Optics code


                  Developed by Michele Carla',
    based on the initial work of Manu Canals and Michele Carla'
           with the support of ALBA (www.cells.es)

                https://github.com/mcarla/ufo


Requirements
------------
The following packages are required to run UFO:

- python 3.5 or later
- numpy
- pyopencl at least one properly configured OpenCL back-end
- ipython3 is required for interactive use only


Install
-------
The latest development version of UFO can be retrieved from github and installed with:

```
git clone https://github.com/mcarla/ufo
pip3 install -e ufo/
```


Getting started
---------------

### Check OpenCL back-ends availability

At least one properly configured OpenCL back-end is required to run any simulation,
a list of the available OpenCL back-ends can be obtained with:
```
import ufo

ufo.list_devices()
```

The output should resemble:

```
0:   Quadro P600
1:   pthread-Intel(R) Core(TM) i5-8400 CPU @ 2.80GHz
```

In this example two back-ends are available: 0 is an Nvidia Quadro GPU, while 1 is an Intel i5 CPU.
UFO can be called from a script or used interactively for example with ipython3. In the latter case some online documentation is accessible via the '?' operator, for example:

```
In [5]: ufo.Quadrupole?
Init signature: ufo.Quadrupole(label, slices=None, length=0.0, k1=0.0, dx=0.0, dy=0.0, dkn=[], dks=[])
Docstring:     
A thick normal quadrupole element.

label    : str   -> Name of the element
slices   : float -> Number of slices used for 'tea pot' expansion (if KICK flag is set)
length   : float -> Length of the element
k1       : float -> Strength of the quadrupolar component
dx, dy   : float -> Horizontal and vertical alignment offset
dkn, dks : list  -> List of normal and skew multipolar field errors
```

How does it work?
-----------------

Besides a few utilities to manipulate the lattice, UFO provides a few classes to run physics simulations (tracking, closed-orbit... ).
Those classes have been designed aiming to perfomance, therefore a few steps required to run the actual simulation
could differ from other simulation codes and may look if anything unusual.
To understand the design choices (and properly use the code) is useful to get an idea of the actual 'behind the scene'.

The entire accelerator lattice, including lines and elements, can be rapresented through the Lattice class.
A lattice object can be build manually or more easily from an input file (UFO uses a MAD-X dialect to describe a lattice).

The lines contained in a Lattice object can then be used to feed one of the simulation classes.
All the physics classes follow the same approach: an OpenCL representation of the `pass method` (routine to track a particle through the line) of the line is produced
and compiled for a specific OpenCL back-end. Subsequently the user can run the simulation as many times as required.

For performance's sake all the elements parameters (e.g. length, gradients) and initial simulation conditions (e.g. initial particle coordinates)
are hardcoded in the OpenCL code allowing the compiler to optimize at best the code.
On the other hand this choice is such that every time the simulation is run the same result is always produced, since there is no way to change any key parameter...

However when creating the simulation the user can provide a list of parameters that will not be hard-coded in the OpenCL code and instead will be specified at run time.
Parameters are specified per particle, therefore it is possible to run in parallel simulations for different 'flavor' of the same lattice,
as required for example to compute a response matrix, where the same simulation is run for different magnets settings.

A few examples
--------------

### Closed orbit

Closed orbit calculation is carried out by the ClosedOrbit class.
In this example the computation is repeated 10 times for 10 different random horizontal misalignment of 'QD'.
note that if 'QD' appears multiple times in the lattice, the same 'dx' error would be applied to the entire family.
To point a specific element, the element-index number can be used (the position in the lattice)

```
import ufo
import numpy

fodo = ufo.Lattice(path='fodo.mad')

count = 10
parameters = [('QD', 'dx')]
orbit = ufo.ClosedOrbit(fodo.RING, particles=count, parameters=parameters)

for i in range(count):
    orbit.parameters[i] = numpy.random.randn() * 1e-3

orbit.run(threads=count)
print(orbit.orbits)
```

### Dynamic Aperture

Dynamic aperture can be simulated through the StableAperture class.
The parameters 'x' and 'y' allow to set the initial particles coordinates.
The output of the simulation is the number of turns each particle survived and can be accessed through the attribute `lost`. 

```
import ufo
import numpy

optics = ufo.Lattice(path='alba.mad')

count = 12 #Simulate a grid of 12x12 particles

parameters = ['x', 'y', ('SF1', 'k2')]
sa = ufo.StableAperture(optics.RING, particles=count**2, turns=1000,
                        flags=ufo.FIVED, parameters=parameters, dp=0.)

x = numpy.linspace(-0.04, 0.04, num=count) #Initial particle coordinates
y = numpy.linspace(-0.04, 0.04, num=count) #Are arranged on a grid
xx, yy = numpy.meshgrid(x, y, sparse=False)

sa.parameters[:, 0] = xx.flatten()
sa.parameters[:, 1] = yy.flatten()
sa.parameters[:, 2] = 25.7971933671 #K2 of SF1

sa.run(threads=count**2)
print('\n        Dynamic aperture with nominal parameters:\n')
print(sa.lost.reshape([count, count]))

sa.parameters[:, 2] = 25.7971933671 * 2.0 #K2 let's double the K2 of SF1
sa.run(threads=count**2)

print('\n        Dynamic aperture with double strength for SF1:\n')
print(sa.lost.reshape([count, count]))
```

### Twiss parameters

In the following example the periodic optics functions for a simple FODO channel are computed.
The positions at which position dependent functions (betatron amplitudes, dispersion...) are evaluated, are specified as a list of float through the keyword argument 'where'. The integer part of each float represent the element index (in this case 1, the second element which happen to be a drift), while the fractional part represents a fraction of the length of the selected elemnt. Therefore in this example the functions will be evaluated in the middle of the first sdrift (to indicate the end of the line the jolly identifier -1 can be used):

```
import ufo

lat = ufo.Lattice(path='fodo.mad')
opt = ufo.Optics(lat.RING, where=[1.5], periodic=True)
opt.run()
```

Results can be accessed through the class attributes: `bx, by, ax, ay, dx, dy, qx, qy` (i.e. `opt.qx`)
Note that when computing periodic solutions optics functions are always computed also at the end of the line.

Supported elements and parameters
---------------------------------

UFO support several elements (drift, quarupole, dipole...). Each element is identified by a unique label,
furthermore the following keyword parameters are supported by each element:


|            | slices | length | angle |  k1   | e1 / e2 | hgap / fint | knl / ksl | k2 / k2s | k3 / k3s | dx / dy | dkn / dks |
|       ---: |  :---: |  :---: | :---: | :---: |  :---:  |    :---:    |   :---:   |  :---:   |  :---:   |  :---:  |   :---:   |
| Marker     |        |        |       |       |         |             |           |          |          |         |           |
| Drift      |        |   x    |       |       |         |             |           |          |          |         |           |
| Mutlipole  |        |        |       |       |         |             |     x     |          |          |    x    |           |
| Quadrupole |    x   |   x    |       |   x   |         |             |           |          |          |    x    |     x     |
| Sbend      |    x   |   x    |   x   |   x   |    x    |      x      |           |          |          |    x    |     x     |
| Rbend      |    x   |   x    |   x   |   x   |    x    |      x      |           |          |          |    x    |     x     |
| Sextupole  |    x   |   x    |       |       |         |             |           |    x     |          |    x    |     x     |
| Octupole   |    x   |   x    |       |       |         |             |           |          |    x     |    x    |     x     |

For Quadrupole, Sbend and Rbend the `slices`, `dkn` and `dks` parameters have effect only when the `KICK` flag is set to True (Teapot expansioin)


**slices:**      Number of slices used in Teapot expansion         

**length:**      Length of the element [m]

**angle:**       Bending angle [rad]

**k1:**          Magnetic gradient normalized by the beam rigidity [m^-2]
 
**e1 / e2:**     Entry / exit pole angle [rad]

**hgap / fint:** Magnet half gap [m] and fringe field integral

**knl / ksl:**   Integrated normal and skew field multipolar components [m^-n]

**k2 / k2s:**    Normal and skew sextupolar component [m^-3]       

**k3 / k3s:**    Normal and skew octupolar component [m^-4]        

**dx / dy:**     Horizontal and vertical alignement errors [m]        

**dkn / dks:** Integrated normal and skew multipolar field errors [m^-n] 

For example a new quadrupole with length 1 m and strength 2 m^-1 can be instantiated with:

```
qf = ufo.Quadrupole('qf', length=1, k1=2)
```

Afterwards, element parameters can be accessed and changed with:

```
qf.dx = 1e-3   #displace horizontally 'qf' by 1mm
qf.k1 *= 1.1   #increase k1 by 10%
```

Field multipoles `knl, ksl, dkn and dks` can be accessed in two diffrent ways

* as a vector:

```
qf.dkn = [0., 0., 0.3]  #Add a sextupolar error to 'qf'
```

* With a scalar notation `kN, kNs, dkN and dkNs` with N the order of the multipole. Note that in this case the component must have been previously defined when the element was instantiated or with the vector notation.
  As explained afterwards, this notation is useful when speciying parameters. Whereas the vector notation is not compatible with parametrization.

```
#if 'knl' has length < 3, the next instruction will fail
qf = ufo.Quadrupole('qf', length=1, k1=2, dkn=[0, 0, 0])

qf.dk2 = 0.3 #has the same effect as qf.dkn = [0., 0., 0.3]
```

Lines
-----

A line is a list of elements or other lines and is used to describe the ordering of a magnetic lattice.
The Line class is derived from the python list class, therefore all the standard lists functionallities are available (`append, remove, copy...`)
The Line class constructor is defined as:
### Line(label, line=[])

**label:**   An unique string identifying the line

**line:**    A list of elements or other lines to be included

A new line can be created with:
```
l = ufo.Line('myline', [qf, qd])
```

The Line class is provided with some special function and attribute to easy certain common tedious tasks:

**count:**  return the number of elements in the line (including the ones in included lines)

**length:**    the overall line length
**angle:**     the overall line bending angle
**flatten():** return a flattened version of the line (all sub-lines will be recursively expanded)
**find():**    find the element index given a filter function provided by the user
**locate():**  return the position in meters of a given element index

