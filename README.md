# real2GAME

The utilities in this repository create initial conditions for [GAME](https://github.com/OpenNWP/GAME) as well as initial and boundary conditions for [L-GAME](https://github.com/OpenNWP/L-GAME).

## Code structure

real2GAME incorporates different tools:

* The interpolator executes the main task: it interpolates the data from another model to the grid of GAME or L-GAME.
* The formatter needs to be executed before the interpolator. It brings the data from the foreign model into a standardized format.
* The downloader needs to be executed before the formatter. It downloads the data from the server of the foreign model.
* The interpolation creator prepares the interpolation indices and weights.

## Building

It is recommended to run real2GAME and (L-)GAME on Linux. These installation instructions are tested for Ubuntu, for other Linux distributions they might have to be modified.

### Dependencies

	sudo apt-get install libeccodes-dev bzip2 bc

Additionally, GAME and/or L-GAME need to be installed, depending on which model you want to use.

### Download and installation

	git clone https://github.com/OpenNWP/real2GAME.git
	cd real2GAME
	./create_directories.sh
	./compile_all.sh -f

## Execution

This is a first try for creating a model input file:

	cd interpolation_creator
	./run.sh
	cd ..
	./run_manually.sh

Modify the directories in `run.sh` and `run_manually.sh` according to your directory structure.
	
