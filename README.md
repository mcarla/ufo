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
- pyopencl with a properly configured OpenCL back-end
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

### Twiss parameters example

In the following example the periodic optics functions for a simple FODO channel are computed.
The positions at which position dependent functions (betatron amplitudes, dispersion...) are evaluated, are specified as a list of float through the keyword argument 'where'. The integer part of each float represent the element index (in this case 1, the second element which happen to be a drift), while the fractional part represents a fraction of the length of the selected elemnt. Therefore in this example the functions will be evaluated in the middle of the first sdrift (to indicate the end of the line the jolly identifier -1 can be used):

```
import ufo

lat = ufo.Lattice(path='fodo.mad')
opt = ufo.Optics(lat.RING, where=[1.5], periodic=True)
opt.run()
opt.qx
```

Results can be accessed through the class attributes: `bx, by, ax, ay, dx, dy, qx, qy`
Note that when computing periodic solutions optics functions are always computed also at the end of the line.

