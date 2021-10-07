
"                            U.F.O."

```
 
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

```

***UFO*** is a fast accelerator optics toolkit designed with GPU in mind, nevertheless it can be used on CPUs with good results.
UFO is not meant to be a general purpose tool, instead it aims to performance at expenses of flexibility and ease of use.

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

Line(label, line=[])
--------------------

A line is a list of elements or other lines and is used to describe the ordering of a magnetic lattice.
The Line class is derived from the python list class, therefore all the standard lists functionalities are available (`append, remove, copy...`)

* **label:**   An unique string identifying the line

* **line:**    A list of elements or other lines to be included

The Line class is provided with some special function and attribute to easy certain common tedious tasks:

Attributes:

* **count:**  Teturn the number of elements in the line (including the ones in included lines)

* **length:** The overall line length

* **angle:**  The overall line bending angle

Methods:

* **flatten():**        return a flattened version of the line (all sub-lines will be recursively expanded)

* **find(what):**       return a list of elements indices such that `what(element) == True`
  * ***what:*** filter function should take as argument an element object and return a boolean (lambda functions here are pretty handy)

* **locate(element):**  return the position in meters of a given element index
  * ***element:*** element object

Examples:
```
import ufo

bend = ufo.Sbend('B', length=1., angle=0.1)
qf = ufo.Quadrupole('QF', length=0.3, k1=1.)
qd = ufo.Quadrupole('QD', length=0.3, k1=-1.)

period = ufo.Line('fodo', [qf, bend, qd, bend])
ring = ufo.Line('ring', [period, period, period, period])

ring.angle
 > 0.8
 
ring.find(lambda e: type(e)==ufo.Quadrupole)
 > [0, 2, 4, 6, 8, 10, 12, 14]

ring.find(lambda e: e.label=='QD')
 > [2, 6, 10, 14]
```

Lattice(path=None)
------------------

A lattice is a collection of elements and lines, it can be assembled manually or more conveniently from a file.
The input file syntax is compatible with the MAD-X syntax, even if not all the MAD-X features are supported.

* **path:** input file path, if not specified an empty lattice is generated

Examples:
```
import ufo

lattice = ufo.Lattice('ring.mad')
lattice.RING.length
 > 10.4
```

where 'ring.mad' is:
```
QF: QUADRUPOLE, L=0.3, K1= 1;
QD: QUADRUPOLE, L=0.3, K1=-1;

B: SBEND, L=1, ANGLE=0.1;

PERIOD: LINE=(QF, B, QD, B);
RING:   LINE=(PERIOD, PERIOD, PERIOD, PERIOD);
```

Beam(energy=3e9, particle_mass=electron_mass, bunch_charge=1e-9, beam_current=0.25, ex=1e-9, ey=1e-9, bunch_length=6e-3, energy_spread=1e-3)
----------------------
A Beam object collects all the informations relative to particles

* **energy:**        particles energy [eV]

* **particle_mass:** particle mass [eV]

* **bunch_charge:**  charge of the entire bunch [C]

* **beam_current:**  average beam current [A]

* **ex / ey:**       horizontal and vertical emittances [m]

* **bunch_length:**  longitudinal dimension of the bunch [m]

* **energy_spread:** RMS energy spread of the bunch divided by the average bunch energy

dump(lattice, path, style='mad', beam=Beam())
---------------------------------------------

A lattice can be exported in a text file compatible with one of the following programs MAD-X, Elegant, Accelerator-Toolbox or OPA

* **lattice:** the Lattice to be exported

* **path:**    the name of the output file to be generated

* **style:**   A string representing the syntax to be used. It can be' mad', 'elegant', 'at' or 'opa'

* **beam:**    Beam informations are used only in when exporting to 'at' or 'opa'

Examples:
```
import ufo

lattice = ufo.Lattice('ring.mad')
ufo.dump(lattice, 'ring.ele', style='elegant')
```

Track(line, flags=0x00, turns=1000, particles=1000, parameters=[], where=[], dp=0., context=None, options=None)
-------------------

Track allows to track a bunch of particles through a line. Tracking parameters includes intial particle conditions and lattice parameters (length of elements, field strenghts, alignment errors...). Parameters are specified of a per-particle basis, therefore it is possible to track at the same time (in parallel hardware permitting) particles with different optics settings, for example with different sextupols strength or different alignment errors.

* **line:**  the Line object to be used for tracking

* **flags:** flags to control specific tracking options, see section **flags** for detailed informations

* **turns:** the number of turns to tracked

* **particles:** the number of particles to be tracked

* **parameters:** a list of parameters to be allowed for variation, see section **parameters** for detailed informations

* **where:** a list of positions where the particles coordinate will be recorded at each turn. Positions are specified as a floating variable, where the integer part identifies the element index while the fractional part is the position along the element expressed as a fraction of the entire element length. The special position -1 identify the end of the line

* **dp:** relative energy deviation used when the **FIVED** flag is set (5D simulation)

* **context:** an OpenCL context as returned by `ufo.context()`

* **options:** OpenCL back-end options, see section **OpenCL options** for detailed informations

Attributes:

* **parameters:** a numpy buffer containing all the simulation parameters. Note that the buffer is not initialized, therefore is up to the user to set it properly before calling the `run()` method

* **tracks** the output buffer where the tracking results will be stored

* **src:** source code of the OpenCL kernel, useful for debugging

Methods:

* **run(threads=1):** run the simulation. The simulation can be run as many times as necessary and the parameters changed between each run
  * **threads:** number of threads to be run in parallel. Is up to the user to determine the best value in terms of performance. For a CPU the optimum is usually the number of available cores. For a GPU the optimum is usually a multiple of the number of cores (2, 3 or 4 times the number of cores seems to be the sweet spot). Therefore for a small GPU threads can be as high as ~10^3, for high end GPU ~10^4 is normal.

Parameters
----------
Parameters are used to identify which variables will be specified independently for each particle and or need to be changed between one simulation run and another.
There is no distinction between beam (intial particle coordinates...) and lattice (field strength...) .
Particle coordinates parameters are identified by the strings: **'x', 'px', 'y', 'py', 'z', 'dp'**. While lattice parameters are identified by tuples of two elements, the first one identify an element or a family of elements and the second one the actual name of the parameter.
When not specified the beam parameters x, px... will be set to 0, while the lattice parameters will be set to the value defined in the lattice at the time the simulation is instantiated. Once a simulation has been instantiated, changing any value of an element will not have any effect on the simulations, only parameters defined throught the parameters argument can be changed (through the parameters attribute). 

Examples:

'k1' of the first element of a line is a parameter:
```
parameters = [(0, 'k1')]
```

'k1' of the element QF is a parameter. If QF appears multiple times (family), every appearence will be affected the same:
```
parameters = [('QF', 'length')]
```

Multiples parameters can be specified for one simulation. In this case we want to define for each particle the initial coordinates x and y, the length of the family QF and the field component k2 of the third element:
```
parameters = ['x', 'y', ('QF', 'length'), (2, 'k2')]
```

Flags
-----

Many simulations parameters are controlled by a set of boolean flags:

* **LINEAR:**           purely linear simulation. Every non linear field will be set to 0. Useful for example for linear optics

* **FIVED:**            particles energy is fixed to a common value and not allowed to vary during the simulation (5D simulation). 5D simulations can be faster respect to full 6D

* **EXACT:**            use exact Hamiltonian

* **KICK:**             replace thick elements (Sbend, Rbend and Quadrupole) with thin kicks and drifts using the teapot expansion. The default number of slices used for the expansion is defined by the variables:
  * **DEFAULT_QUADRUPOLE_SLICES**
  * **DEFAULT_BEND_SLICES**
  * **DEFAULT_MULTIPOLE_SLICES**
  * **DEFAULT_SEXTUPOLE_SLICES**
  * **DEFAULT_OCTUPOLE_SLICES**

* **RADIATION:**        not yet implemented

* **DOUBLE_PRECISION:** use 64 bit variables instead of 32 bit. Double precision allows for higher precision at the expense of performance (expecially on GPU)

* **ACHROMATIC:**       suppress the chromatic focusing effects. Useful for dispersion computation

OpenCL options
--------------
By default the following OpenCL options are passed to the OpenCL back-end:

* -cl-fast-relaxed-math
* -cl-mad-enable
* -cl-single-precision-constant

The latter is removed when switching to 64 bit variables (when the flag **DOUBLE_PRECISION** is set). The default options can be changed by setting the string: `ufo.DEFAULT_CL_OPTIONS`
