# real2GAME

The utilities in this repository create initial conditions for [GAME](https://github.com/OpenNWP/GAME) as well as initial and boundary conitions for [L-GAME](https://github.com/OpenNWP/L-GAME).

## Building

### Dependencies

* [geos95](https://github.com/OpenNWP/geos95)
* netcdf library (Ubuntu: `sudo apt-get install libnetcdf-dev`)
* bzip2 (Ubuntu: `sudo apt-get install bzip2`)
* basic calculator (Ubuntu: `sudo apt-get install bc`)

#### For developing

* Valgrind (Ubuntu: `sudo apt-get install valgrind`, for doing checks)

### Download and installation

	git clone https://github.com/OpenNWP/real2GAME.git
	cd real2GAME
	./compile_all.sh
