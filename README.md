
```
                   ██    ██ ███████  ██████  
                   ██    ██ ██      ██    ██ 
                   ██    ██ █████   ██    ██ 
                   ██    ██ ██      ██    ██ 
                    ██████  ██       ██████ 

 
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

...***UFO*** is a fast accelerator optics toolkit designed with GPU in mind, nevertheless it gets along well with CPUs too.
UFO is not meant to be a general purpose tool, instead it aims to performance at expenses of flexibility and ease of use...

## Requirements

The following packages are required to run UFO:

- python 3.5 or later
- numpy
- pyopencl at least one properly configured OpenCL back-end
- ipython3 is required for interactive use only


## Install

The latest development version of UFO can be retrieved from github and installed with:

```
git clone https://github.com/mcarla/ufo
pip3 install -e ufo/
```


## Getting started

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

Refer to the [HOWTO](https://github.com/mcarla/ufo/blob/main/HOWTO.md) for further documentation.
