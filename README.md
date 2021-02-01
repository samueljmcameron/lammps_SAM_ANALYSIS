# lammps_SAM_ANALYSIS

Custom LAMMPS code for analysing LAMMPS data. To use, copy and paste this
repo into the LAMMPS source code (which can be downloaded from e.g.
https://github.com/lammps/lammps). This repo is up-to-date with LAMMPS
stable release from 29 October 2020. Once copied and pasted, add the package
to the cmake list optional packages by modifying the line in
/path/to/lammps/cmake/CMakeLists.txt which sets all the packages list
(should look like "set(STANDARD_PACKAGES ASPHERE ... USER-ADIOS)" for stable
release 29 October 2020, and starts on line 107) to include SAM-ANALYSIS
(so the above line should now look like "set(STANDARD_PACKAGES ASPHERE ...
USER-ADIOS SAM-ANALYSIS)" for stable release 29 October 2020).

Then, make lammps using cmake, with the flag PKG_SAM-ANALYSIS=yes set:

cmake -D PKG_SAM-ANALYSIS=yes /path/to/lammps/cmake
make

and you'll have all the functionality of this repo integrated into lammps.
Note that a substantial portion of this repo uses the dipole package, so
you'll probably want to add the PKG_DIPOLE=yes flag to cmake as well.

Once you've done the above, you should have all the functionality of
LAMMPS plus these additional files. The documentation for the
custom LAMMPS computes in this repo should also show up if you make the
docs, though you'll have to search for each html page individually instead
of just e.g. looking at the list of computes in
/path/to/lammps/doc/html/Commands_compute.html, since fully integrating
the custom docs is a bit of a headache.


One warning: When doing the above, this repo will overwrite the files
src/read_dump.cpp, and src/reader_native.cpp with exact duplicates (for
stable release 29 October 2020) aside from adding a few lines of code to
allow for reading in dipole data from dump files. This is necessary for
the custom compute rdf/dipole if one wants to read data from dump files
using the rerun command of lammps.