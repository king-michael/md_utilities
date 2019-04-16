# VMD utilities

## connection_python
Establish a connection between python
and VMD via sending data over a socket.

## gopython

### live_fes.py
#### usage:
First load the script, then use ``run()`` to load and connect your files.

**Load the script**:
```
gopython /home/micha/SIM-PhD-King/gits/md_utilites/vmd_util/live_fes.py
```

**start the command**:
```
gopython -command "run('analysis/fes.dat', 'analysis/COLVAR')"
```

#### important functions

* `run(fes_file='fes.dat', colvar='COLVAR')` <br>
    Plots a `fes.dat` file and
    connects the current ``vmd_frame`` with the ``COLVAR`` file.
    This allows to live trace where the structe currently is.