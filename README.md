3D and 4D variational data assimilation for a wide range of geophysical models.

# Building

## Dependencies

* [geos95](https://github.com/MHBalsmeier/geos95)
* netcdf library (Ubuntu: sudo apt-get libnetcdf-dev)
* CMake (Ubuntu: sudo apt-get install cmake)
* [atmostracers](https://github.com/MHBalsmeier/atmostracers)
* OpenMPI (Ubuntu: sudo apt-get install mpich)

### For developing

* Valgrind (Ubuntu: sudo apt-get install valgrind, for doing checks)

## Download

	git clone https://github.com/MHBalsmeier/ndvar.git
	cd ndvar

## Compilation

In order to compile the model, cd into the directory build\_scripts and execute the shell script build\_install.sh. Scripts with the suffix \_dev are installing to a destination where new versions can be tested.
