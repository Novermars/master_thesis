# Cerebral aneurysm modelling using the Lattice Boltzmann method

* This is the main repository for my master thesis, on different possible collision operators in the lattice Boltzmann method using the LBM package waLBerla.

## Build instructions:
1) Make a new directory build and cd into that directory
2) `cmake ..` (You may have to manually specifiy the location of waLBerla on your system)
3) Be sure to turn on code_gen and building with openmesh via e.g. ccmake, turn off tutorials/benchmarks/showcases
4) make -j number_of_cores
5) ./src/aneurysm_simulation src/config.prm

You can use Paraview to look at the results.
