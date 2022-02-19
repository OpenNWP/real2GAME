# real2GAME

The utilities in this repository create initial conditions for [GAME](https://github.com/OpenNWP/GAME) as well as initial and boundary conitions for [L-GAME](https://github.com/OpenNWP/L-GAME).

## Code structure

real2GAME incorporates different tool

* The interpolator executes the main task: it interpolates the data from another model to the grid of GAME or L-GAME.
* The formatter needs to be executed before the interpolator. It brings the data from the foreign model into a standardized format.
* The downloader needs to be executed before the formatter. It downloads the data from the server of the foreign model.
* The interpolation creator prepares the interpolation indices and weights to increase the performance.

## Building

### Dependencies

* [GAME](https://github.com/OpenNWP/GAME)
* bzip2 (Ubuntu: `sudo apt-get install bzip2`)
* basic calculator (Ubuntu: `sudo apt-get install bc`)

### Download and installation

	git clone https://github.com/OpenNWP/real2GAME.git
	cd real2GAME
	./create_directories.sh
	./compile.sh
