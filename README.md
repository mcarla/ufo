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
- pyopencl with a properly configured OpenCL backend
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

At least one properly configured OpenCL backend is required to run any simulation,
a list of the available OpenCL backends can be obtained with:
```
import ufo

ufo.list_devices()
```

The output should resemble:

```
0:   Quadro P600
1:   pthread-Intel(R) Core(TM) i5-8400 CPU @ 2.80GHz
```

In this case two backends are available: 0 is an envidia Quadro GPU, while 1 is an Intel i5 CPU.

