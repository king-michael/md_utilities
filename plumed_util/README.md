# PLUMED utilities


## switchingfunction
Various switching functions used in ``PLUMED`` as described under
[switchingfunction](https://plumed.github.io/doc-master/user-doc/html/switchingfunction.html)

* ``rational(r, r_0, d_0=0.0, n=6, m=None, d_max=None)`` <br>
    RATIONAL switching function from PLUMED

    `s(r) = ( 1 - ( (r-d_0) / r_0 )^n ) / ( 1 - ( (r-d_0) / r_0 )^m )`


## fes_file
Functions to load a ``fes.dat`` file produced by ``plumed sum_hills``.

* ``get_header(inp_file)`` <br>
    Returns the header of a `fes.dat` file.

* ``load(inp_file)`` <br>
    Function to load the ``fes.dat`` file produced by ``PLUMED``.

* ``plot(header, data, ax=None)`` <br>
    Function to plot a ``fes.dat`` file.