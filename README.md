nxsunit
=======

Data field loader and unit conversion for NeXus HDF files.

NeXus data files from various institutions can use a wide variety of units.
This package attempts to capture the various unit definitions and provide
conversions so that the user can load data in consistent units. For example,
all instrument geometry can be loaded in meters even if some components are
stored in centimeters and others in millimeters.


See `nxsunit.py` for interface details.


Change history
==============

1.0 2022-07-12
----------------

* Unify units code from scattering/scattio, sasview, reductus, and nexusformat
* Add the `data_as` function from reductus for loading data fields
* Create pip installable package
