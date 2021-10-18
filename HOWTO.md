
## _How does it work?_

Besides a few utilities to manipulate the lattice, UFO provides a few classes to run physics simulations (tracking, closed-orbit... ).
Those classes have been designed aiming to perfomance, therefore a few steps required to run the actual simulation
could differ from other simulation codes and may look if anything unusual.
To understand the design choices (and properly use the code) is useful to get an idea of the actual 'behind the scene'.

The entire accelerator lattice, including lines and elements, can be rapresented through the Lattice class.
A lattice object can be build manually or more easily from an input file (UFO uses a MAD-X dialect to describe a lattice).

The lines contained in a Lattice object can then be used to feed one of the simulation classes.
All the physics classes follow the same approach: an OpenCL representation of the 'pass method' (routine to track a particle through the line) of the line is produced
and compiled for a specific OpenCL back-end. Subsequently the user can run the simulation as many times as required.

For performance's sake all the elements parameters (e.g. length, gradients) and initial simulation conditions (e.g. initial particle coordinates)
are hardcoded in the OpenCL code allowing the compiler to optimize at best the code.
On the other hand this choice is such that every time the simulation is run the same result is always produced, since there is no way to change any key parameter...

However when creating the simulation the user can provide a list of parameters that will not be hard-coded in the OpenCL code and instead will be specified at run time.
Parameters are specified per particle, therefore it is possible to run in parallel simulations for different 'flavor' of the same lattice,
as required for example to compute a response matrix, where the same simulation is run for different magnets settings.


## _Elements_

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


* **slices:**      Number of slices used in Teapot expansion
* **length:**      Length of the element [m]
* **angle:**       Bending angle [rad]
* **k1:**          Magnetic gradient normalized by the beam rigidity [m^-2]
* **e1 / e2:**     Entry / exit pole angle [rad]
* **hgap / fint:** Magnet half gap [m] and fringe field integral
* **knl / ksl:**   Integrated normal and skew field multipolar components [m^-n]
* **k2 / k2s:**    Normal and skew sextupolar component [m^-3]       
* **k3 / k3s:**    Normal and skew octupolar component [m^-4]        
* **dx / dy:**     Horizontal and vertical alignement errors [m]        
* **dkn / dks:** Integrated normal and skew multipolar field errors [m^-n] 

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

New elements can be defined from scratch or derived from existing one. A typical use case is defining elements with time dependent behaviour (e.g. kicker, power supply noise...). In this case the only part of the element class that needs to be modified is the 'pass method' which is called `method()`.
Follows a simple simulation of the emittance blow-up induced by a quadrupolar kick, as the one produced on the stored beam by the residual field of a non-linear injection kicker.
The non-linear kicker is simulated with a thin multipolar multipole, the 'pass method' is redefined in order to switch on the kick only at turn 20:

```
import numpy
import ufo

class NLK(ufo.Multipole):
    def method(self, k1=None, **kwargs):
        k1 = k1 or self.knl[1]
        k1 = f"(turn == 20) ? {k1} : 0."

        return super().method(k1=f'({k1})', **kwargs)

alba = ufo.Lattice(path='../optics/alba.mad')
nlk = NLK('NLK', knl=[0., 5e-2])

alba.RING.insert(0, nlk)

betax = 11.26
ex0 = 4.3e-9

count = 100
track = ufo.Track(alba.RING, particles=count, flags=ufo.FIVED, turns=65, where=[-1], parameters=['x', 'px'])

numpy.random.seed(0)
track.parameters[:, 0] = numpy.random.randn(count) * numpy.sqrt(0.5 * betax * ex0)
track.parameters[:, 1] = numpy.random.randn(count) * numpy.sqrt(0.5 * ex0 / betax)

track.run(threads=100)
```

Results are plotted with:

```
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
ex = track.tracks[:, :, 0]**2 / betax + track.tracks[:, :, 1]**2 * betax
ex = numpy.mean(ex, axis=0)
ax.plot(ex * 1e9)
ax.set_xlabel('Turn #')
ax.set_ylabel('Ex [nm]')
plt.show()
```

![nlk emittance blow-up](https://github.com/mcarla/ufo/blob/main/doc/nlk.png)

## _Line(label, line=[])_

A line is a list of elements or other lines and is used to describe the ordering of a magnetic lattice.
The Line class is derived from the python list class, therefore all the standard lists functionalities are available (`append, remove, copy...`)

* **label:**   An unique string identifying the line
* **line:**    A list of elements or other lines to be included

The Line class is provided with some special function and attribute to easy certain common tedious tasks:

### Attributes:

* **count:**  Teturn the number of elements in the line (including the ones in included lines)
* **length:** The overall line length
* **angle:**  The overall line bending angle

### Methods:

* **flatten():**        return a flattened version of the line (all sub-lines will be recursively expanded)
* **find(what):**       return a list of elements indices that satisfy the condition `what(element) == True`, where what() is a user specified test function that takes as argument an element object and return a boolean (lambda functions here are pretty handy)
* **locate(element):**  return the position in meters of a given element index. The element index is specified as a floating variable, where the integer part identifies the element index while the fractional part is the position along the element expressed as a fraction of the entire element length. The special position -1 identify the end of the line

### Examples:
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

## _Lattice(path=None)_

A lattice is a collection of elements and lines, it can be assembled manually or more conveniently from a file.
The input file syntax is compatible with the MAD-X syntax, even if not all the MAD-X features are supported.

* **path:** input file path, if not specified an empty lattice is generated

### Examples:
```
import ufo

lattice = ufo.Lattice('ring.mad')
lattice.RING.length
 > 10.4
```

where the input file 'ring.mad' is:
```
QF: QUADRUPOLE, L=0.3, K1= 1;
QD: QUADRUPOLE, L=0.3, K1=-1;

B: SBEND, L=1, ANGLE=0.1;

PERIOD: LINE=(QF, B, QD, B);
RING:   LINE=(PERIOD, PERIOD, PERIOD, PERIOD);
```

## _Beam(energy=3e9, particle_mass=electron_mass, bunch_charge=1e-9, beam_current=0.25, ex=1e-9, ey=1e-9, bunch_length=6e-3, energy_spread=1e-3)_

A Beam object collects all the informations relative to particles

* **energy:**        particles energy [eV]
* **particle_mass:** particle mass [eV]
* **bunch_charge:**  charge of the entire bunch [C]
* **beam_current:**  average beam current [A]
* **ex / ey:**       horizontal and vertical emittances [m]
* **bunch_length:**  longitudinal dimension of the bunch [m]
* **energy_spread:** RMS energy spread of the bunch divided by the average bunch energy

## _dump(lattice, path, style='mad', beam=Beam())_

A lattice can be exported in a text file compatible with one of the following programs MAD-X, Elegant, Accelerator-Toolbox or OPA

* **lattice:** the Lattice to be exported
* **path:**    the name of the output file to be generated
* **style:**   A string representing the syntax to be used. It can be 'mad', 'elegant', 'at' or 'opa'
* **beam:**    Beam informations are used only in when exporting to 'at' or 'opa'

### Examples:
```
import ufo

lattice = ufo.Lattice('ring.mad')
ufo.dump(lattice, 'ring.ele', style='elegant')
```

## _StableAperture(line, flags=0x00, turns=1000, particles=1000, parameters=[], dp=0., context=None, options=None)_

Compute how many turns a particle survive in a ring. At the end of each turn the horizontal and vertical coordinates are checked, if found greater than 1m the particle will be considered as lost. Useful to estimate dynamic aperture and momentum aperture. Simulation parameters includes initial particle coordinates and lattice parameters (length of elements, field strenghts, alignment errors...). Parameters are specified per-particle, therefore it is possible to track at the same time (in parallel hardware permitting) particles with different optics settings, for example with different sextupols strength or different alignment errors.

* **line:**  the Line object to be used for tracking
* **flags:** flags to control specific tracking options, see section **flags** for detailed informations
* **turns:** the maximum number of turns to be tracked
* **particles:** the number of particles to be tracked
* **parameters:** a list of parameters to be allowed for variation, see section **parameters** for detailed informations
* **dp:** relative energy deviation used when the **FIVED** flag is set (5D simulation)
* **context:** an OpenCL context as returned by `ufo.context()`
* **options:** OpenCL back-end options, see section **OpenCL options** for detailed informations

### Attributes:

* **parameters:** a numpy buffer containing all the simulation parameters. Note that the buffer is not initialized, therefore is up to the user to set it properly before calling the `run()` method
* **lost** the output buffer where the number of turns made by each particle before getting lost will be stored
* **src:** source code of the OpenCL kernel, useful for debugging

### Methods:

* **run(threads=1):** run the simulation. The simulation can be run as many times as necessary and the parameters changed between each run
  * **threads:** number of threads to be run in parallel. Is up to the user to determine the best value in terms of performance. For a CPU the optimum is usually the number of available cores. For a GPU the optimum is usually a multiple of the number of cores (2, 3 or 4 times the number of cores seems to be the sweet spot). Therefore for a small GPU threads can be as high as ~10^3, for high end GPU ~10^4 is normal.

### Examples:

1. Compute the dynamic aperture of the ALBA ring by tracking a grid of 12x12 particles equally spaced in the transverse plane:
```
import numpy
import ufo

optics = ufo.Lattice(path='../optics/alba.mad')

count = 12

parameters = ['x', 'y'] #the only simulation parameters are the horizontal and vertical inital position of each particle
sa = ufo.StableAperture(optics.RING, particles=count**2, turns=1000, flags=ufo.FIVED, parameters=parameters)

x = y = numpy.linspace(-0.04, 0.04, num=count)
xx, yy = numpy.meshgrid(x, y, sparse=False)

sa.parameters[:, 0] = xx.flatten()   #initialize the horizontal and
sa.parameters[:, 1] = yy.flatten()   #the vertical initial position of each particle

sa.run(threads=count**2)  #start the simulation running 144 threads in parallel
print(sa.lost.reshape([count, count]))  #print how many turns each particle survived
```
The attribute 'lost' contains the number of turns each particle survived in the ring (in this particular case 1000 is the maximum number of tracked turns).
This examample code should produce the following output:


```
[[   0    0    0    0    0    0    0    0    0    0    0    0]
 [   0    0    0    0    0    0    0    0    0    0    0    0]
 [   0    1    1    2    0    0    1    1    2    2    1    0]
 [   0    2   30    4    8    1    7   14    8    2    6    0]
 [   0    1   11 1000 1000 1000 1000  981   52    6    1    0]
 [   0    4 1000 1000 1000 1000 1000 1000 1000 1000    2    1]
 [   0    4 1000 1000 1000 1000 1000 1000 1000 1000    2    1]
 [   0    1   11 1000 1000 1000 1000  981   52    6    1    0]
 [   0    2   30    4    8    1    7   14    8    2    6    0]
 [   0    1    1    2    0    0    1    1    2    2    1    0]
 [   0    0    0    0    0    0    0    0    0    0    0    0]
 [   0    0    0    0    0    0    0    0    0    0    0    0]]
```

2. An octupolar field error is now added independently to each quadrupole.
Note that field errors (dkn / dks) are not implemented in thick elements such as Sbend, Rbend and Quadrupole, therefore the teapot expansion should be enabled by setting the **KICK** flag when calling StableAperture. The teapot expansion will convert thick elements in a succession of thin kicks and drifts, allowing to include higher order multipoles but resulting in a usually slower execution time. More informations on how to set the number of slices used in the teapot expansion is found in the section **flags**.
```
import numpy
import ufo

alba = ufo.Lattice(path='alba.mad')

count = 12 #will be tracking a grid of 12x12 particles

#we want to set an octupolar error at each quadrupole -> find each individual quadrupole:
quadrupoles  = alba.RING.find(lambda e: type(e) == ufo.Quadrupole)

parameters  = ['x', 'y'] #iniatl particle coordinates
parameters += [(e, 'dk3') for e in quadrupoles] #add individual octupolar errors at each quadrupole

sa = ufo.StableAperture(alba.RING, particles=count**2, turns=1000, flags=ufo.FIVED | ufo.KICK, parameters=parameters)

x = y = numpy.linspace(-0.04, 0.04, num=count) #initial particle coordinates
xx, yy = numpy.meshgrid(x, y, sparse=False)    #are on a square grid 4x4 cm^2

sa.parameters[:, 0] = xx.flatten() #set the inial particles coordinates
sa.parameters[:, 1] = yy.flatten()

#set a random octupolar error for each magnet, but equal for each particle
numpy.random.seed(0)
sa.parameters[:, 2:] = numpy.random.randn(len(parameters) - 2) * 100.

sa.run(threads=count**2)
print(sa.lost.reshape([count, count]))
```


## _Track(line, flags=0x00, turns=1000, particles=1000, parameters=[], where=[], dp=0., context=None, options=None)_

Track allows to track a bunch of particles through a line. Simulation parameters includes initial particle coordinates and lattice parameters (length of elements, field strenghts, alignment errors...). Parameters are specified per-particle, therefore it is possible to track at the same time (in parallel hardware permitting) particles with different optics settings, for example with different sextupols strength or different alignment errors.

* **line:**  the Line object to be used for tracking
* **flags:** flags to control specific tracking options, see section **flags** for detailed informations
* **turns:** the number of turns to be tracked
* **particles:** the number of particles to be tracked
* **parameters:** a list of parameters to be allowed for variation, see section **parameters** for detailed informations
* **where:** a list of positions where the particles coordinate will be recorded at each turn. Positions are specified as a floating variable, where the integer part identifies the element index while the fractional part is the position along the element expressed as a fraction of the entire element length. The special position -1 identify the end of the line
* **dp:** relative energy deviation used when the **FIVED** flag is set (5D simulation)
* **context:** an OpenCL context as returned by `ufo.context()`
* **options:** OpenCL back-end options, see section **OpenCL options** for detailed informations

### Attributes:

* **parameters:** a numpy buffer containing all the simulation parameters. Note that the buffer is not initialized, therefore is up to the user to set it properly before calling the `run()` method
* **tracks** the output buffer where the tracking results will be stored
* **src:** source code of the OpenCL kernel, useful for debugging

### Methods:

* **run(threads=1):** run the simulation. The simulation can be run as many times as necessary and the parameters changed between each run
  * **threads:** number of threads to be run in parallel. Is up to the user to determine the best value in terms of performance. For a CPU the optimum is usually the number of available cores. For a GPU the optimum is usually a multiple of the number of cores (2, 3 or 4 times the number of cores seems to be the sweet spot). Therefore for a small GPU threads can be as high as ~10^3, for high end GPU ~10^4 is normal.

## _ClosedOrbit(line, flags=0x00, particles=1000, parameters=[], dp=0., iterations=200, context=None, options=None)_

Find the closed orbit by minimizing the residuals between initial and final coordinate of a particle tracked over one turn. The minimization is carried out using the Nelder–Mead method. The optimization process is iterated for a fixed number of times before stopping, no other stopping criteria is available. Simulation parameters includes initial closed orbit 'guess' and lattice parameters (length of elements, field strenghts, alignment errors...). Parameters are specified per-particle, therefore it is possible to compute at the same time (in parallel hardware permitting) closed orbits with different optics settings, useful for example for response matrix.

* **line:**  the Line object to be used for closed orbit computation
* **flags:** flags to control specific tracking options, see section **flags** for detailed informations
* **particles:** the number of closed orbits to be computed
* **parameters:** a list of parameters to be allowed for variation, see section **parameters** for detailed informations
* **dp:** relative energy deviation used when the **FIVED** flag is set (5D simulation)
*  **iteration:** the number of iterations of the Nelder-Mead minimization before stopping the search 
* **context:** an OpenCL context as returned by `ufo.context()`
* **options:** OpenCL back-end options, see section **OpenCL options** for detailed informations

### Attributes:

* **parameters:** a numpy buffer containing all the simulation parameters. Note that the buffer is not initialized, therefore is up to the user to set it properly before calling the `run()` method
* **orbits** the output buffer where the closed orbits will be stored
* **src:** source code of the OpenCL kernel, useful for debugging

### Methods:

* **run(threads=1):** run the simulation. The simulation can be run as many times as necessary and the parameters changed between each run
  * **threads:** number of threads to be run in parallel. Is up to the user to determine the best value in terms of performance. For a CPU the optimum is usually the number of available cores. For a GPU the optimum is usually a multiple of the number of cores (2, 3 or 4 times the number of cores seems to be the sweet spot). Therefore for a small GPU threads can be as high as ~10^3, for high end GPU ~10^4 is normal.

### Examples:

Compute the orbit distortion induced by each horizontal orbit corrector (horizontal orbit response matrix) in the ALBA storage ring:
```
import numpy
import ufo

alba = ufo.Lattice(path='alba.mad')

#ALBA has orbit correctors integrated in the sextupoles
correctors  = alba.RING.find(lambda e: type(e) == ufo.Sextupole)
count = len(correctors) #number of correctors
parameters = [(e, 'dx') for e in correctors] #add horizontal orbit correctors

co = ufo.ClosedOrbit(alba.RING, particles=count, flags=ufo.FIVED, parameters=parameters)
#the computation of the closed orbit distortion due to a corrector is assigned to a particle
co.parameters[:, :] = numpy.identity(count) * 1e-3 
co.run(threads=count) #the computation is run in parallel

print(co.orbits) #print the resulting closed orbit distortion associated to each corrector
```
co.orbits contains the computed closed orbits coordinates (x, px, y, py) at the start of the ring. Each line refers to a different corrector.
The output should look like:
```
[[-3.5447621e-05  2.8344164e-07  3.0312110e-11  1.0759697e-11]
 [ 3.4938712e-05 -3.2280422e-08 -9.1950558e-11 -3.7141124e-11]
 [-2.1214004e-05  1.2076521e-06  1.7580287e-10  9.8768577e-11]
 [ 7.0103335e-05 -2.6784724e-06  1.3894609e-11  3.1442112e-12]
 [-2.8944047e-05  6.6649659e-07 -1.3028309e-10 -3.6299064e-11]
 [ 2.6881035e-05 -3.5004368e-06  1.3336779e-10  3.4140406e-11]
 [-4.7772559e-05  4.6686437e-06  5.3740880e-11 -1.1813918e-12]
 [-4.0694591e-05 -5.1691532e-06 -1.2420462e-10 -2.4220487e-11]
 .....
```
Computation precision can be increased at expenses of speed, by setting the **DOUBLE_PRECISION** flag:
```
co = ufo.ClosedOrbit(alba.RING, particles=count, flags=ufo.FIVED | ufo.DOUBLE_PRECISION, parameters=parameters)
```
in this case the output would be:
```
[[-3.54477825e-05  2.83433330e-07  2.25021296e-14  9.49363972e-15]
 [ 3.49391547e-05 -3.22556427e-08  3.01302730e-14  9.52440403e-15]
 [-2.12181378e-05  1.20743209e-06  1.33593614e-15 -1.50655993e-16]
 [ 7.01036233e-05 -2.67849563e-06  1.72146152e-15  2.17456827e-15]
 [-2.89429422e-05  6.66562989e-07  3.49636411e-14  2.38165648e-15]
 [ 2.68826620e-05 -3.50035798e-06  9.86471226e-15  2.30303455e-15]
 [-4.77740103e-05  4.66845109e-06 -5.43795223e-14 -1.48958729e-14]
 [-4.06958218e-05 -5.16921708e-06  2.87191971e-15 -1.74380874e-14]
 .....
```
Note how in the latter case the vertical coordinates are much closer to 0, as to be expected since no vertical orbit is produced by the correctors and no other errors are present in the lattice.

## _Optics(line, where=[], periodic=True, parameters=[], flags=LINEAR | ACHROMATIC, context=None, options=None_

### Examples:
The positions at which position dependent functions (betatron amplitudes, dispersion...) are evaluated, are specified as a list of float through the keyword argument 'where'. The integer part of each float represent the element index (in this case 1, the second element which happen to be a drift), while the fractional part represents a fraction of the length of the selected elemnt. Therefore in this example the functions will be evaluated in the middle of the first sdrift (to indicate the end of the line the jolly identifier -1 can be used):

```
import ufo

lat = ufo.Lattice(path='fodo.mad')
opt = ufo.Optics(lat.RING, where=[1.5], periodic=True)
opt.run()
```

Results can be accessed through the class attributes: `bx, by, ax, ay, dx, dy, qx, qy` (i.e. `opt.qx`)
Note that when computing periodic solutions optics functions are always computed also at the end of the line.

## _Parameters_

Parameters are used to identify which variables will be specified independently for each particle and or need to be changed between one simulation run and another.
There is no distinction between beam (intial particle coordinates...) and lattice (field strength...) .
Particle coordinates parameters are identified by the strings: `'x', 'px', 'y', 'py', 'z', 'dp'`. While lattice parameters are identified by tuples of two elements, the first one identify an element or a family of elements and the second one the actual name of the parameter.
When not specified the beam parameters x, px... will be set to 0, while the lattice parameters will be set to the value defined in the lattice at the time the simulation is instantiated. Once a simulation has been instantiated, changing any value of an element will not have any effect on the simulations, only parameters defined throught the parameters argument can be changed (through the parameters attribute). 

### Examples:

'k1' of the first element of the line is a parameter:
```
parameters = [(0, 'k1')]
```

'k1' of the element QF is a parameter. If QF appears multiple times (family) in the line, every appearence will be affected the same:
```
parameters = [('QF', 'length')]
```

Multiples parameters can be specified for one simulation. In this case we want to define for each particle the initial coordinates x and y, the length of the family QF and the field component k2 of the third element in the line:
```
parameters = ['x', 'y', ('QF', 'length'), (2, 'k2')]
```

The parameters buffer (simulation attribute) contains one entry for each specified parameter reflecting the ordering of the parameters list,
also the results returned by the simulation (tracking, closed-orbit...) are organized the same way.

### Examples:
```
parameters = [('QD', 'dx'), ('QF', 'dy')]

co = ufo.ClosedOrbit(lattice.RING, particles=3, parameters=parameters) #Simulate 3 particles

co.parameters[0, 0] = 1e-4 #Set QD.dx = 100um for particle 0
co.parameters[1, 0] = 2e-4 #Set QD.dx = 200um for particle 1
co.parameters[2, 0] = 3e-4 #Set QD.dx = 300um for particle 2

co.parameters[:, 1] = 5e-4 #Set QF.dy = 500um for all particles

co.run(threads=10)

print(co.orbits[0]) #Print the orbit of the first particle 
print(co.orbits[1]) #Print the orbit of the second particle
```

## _Flags_

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

## list_devices(device=None)
Print on the screen a list of the available OpenCL back-ends and relative id

* **device** id of a device, if specified will return a device descriptor

#### Examples:
```
import ufo

ufo.list_devices()
```

The output should resemble:
```
0:   Quadro P600
1:   pthread-Intel(R) Core(TM) i5-8400 CPU @ 2.80GHz
```

In this case two back-ends are available. The first device (with id=0) is an Nvidia Quadro P600 GPU, while the second one is an Intel i5 CPU

## _context(device)_
Given a device id it returns an OpenCL context suitable for all the simulation classes (Track, ClosedOrbit...). A list if the available back-ends and relative ids can be obtained with `list_devices()`

* **device** is the OpenCL back-end id (where the first back-end has id 0, the second 1...).

### Examples:
```
import ufo

ctx = ufo.context(0)

lattice = ufo.Lattice(path='alba.mad')

orbit = ufo.ClosedOrbit(lattice.RING, particles=count, context=ctx)
```

## _OpenCL options_

By default the following OpenCL options are passed to the OpenCL back-end:

* -cl-fast-relaxed-math
* -cl-mad-enable
* -cl-single-precision-constant

The latter is removed when switching to 64 bit variables (when the flag **DOUBLE_PRECISION** is set). The default options can be changed by setting the string: `ufo.DEFAULT_CL_OPTIONS`
