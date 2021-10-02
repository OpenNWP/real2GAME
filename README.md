# GAME-DA

Data assimilation for [GAME](https://github.com/OpenNWP/GAME).

## Documents

### Scientific derivation

The scientific derivations of the assimilation techniques can be found in my textbook on theoretical meteorology (in German): [Kompendium Theoretische Meteorologie](https://raw.githubusercontent.com/MHBalsmeier/kompendium/master/kompendium.pdf).

### Documentation

The documentation of the code can be found in the subdirectory doc.

## Building

### Dependencies

* [geos95](https://github.com/OpenNWP/geos95)
* netcdf library (Ubuntu: `sudo apt-get install libnetcdf-dev`)
* [atmostracers](https://github.com/OpenNWP/atmostracers)
* bzip2 (Ubuntu: `sudo apt-get install bzip2`)

#### For developing

* Valgrind (Ubuntu: `sudo apt-get install valgrind`, for doing checks)

### Download and installation

	git clone https://github.com/OpenNWP/GAME-DA.git
	cd GAME-DA
	./build.sh
