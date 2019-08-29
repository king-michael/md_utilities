# lammps_utils
various utilities around [LAMMPS](https://lammps.sandia.gov/)

## Functions

### `conversions.py`
* `convert_box(...)`
    Function to convert box formats

## `units`
* `constants` 
    Dictionary `{ unit : { constant : value } }` for the different 
    [unit](https://lammps.sandia.gov/doc/units.html) used in LAMMPS

### `utils.py`
* `extract_unit_constants_from_src(...)` 
    Function to extract the constants used in LAMMPS formulas from the LAMMPS source files
