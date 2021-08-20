Optimum interpolation (OI), 3D and 4D variational data assimilation for [GAME](https://github.com/OpenNWP/game).

# Documents

## Scientific derivation

The scientific derivations of the assimilation techniques can be found in my textbook on theoretical meteorology (in German): [Kompendium Theoretische Meteorologie](https://raw.githubusercontent.com/MHBalsmeier/kompendium/master/kompendium.pdf).

## Documentation

The documentation of the code can be found in the subdirectory doc.

# Building

## Dependencies

* [geos95](https://github.com/OpenNWP/geos95)
* netcdf library (Ubuntu: `sudo apt-get libnetcdf-dev`)
* CMake (Ubuntu: `sudo apt-get install cmake`)
* [atmostracers](https://github.com/OpenNWP/atmostracers)
* MPICH (Ubuntu: `sudo apt-get install mpich`)
* bzip2 (Ubuntu: `sudo apt-get install bzip2`)
* python3 (Ubuntu: `sudo apt-get install python3`)

### For developing

* Valgrind (Ubuntu: `sudo apt-get install valgrind`, for doing checks)

## Download and installation

	git clone https://github.com/OpenNWP/ndvar.git
	./compile.sh
	cd formatter
	./compile.sh
