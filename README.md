# Cerebral aneurysm modelling using the Lattice Boltzmann method

* This is the main repository for my master thesis, on using the lattice Boltzmann method to simulate bloodflow using the LBM package waLBerla.

* Due to a missing feature in waLBerla's code generation, one of the generated files has to be manually adjusted when newly generated. To circumvent this, this file is already in the git repo. After running the instructions below, you can simply use `git checkout .` and the proper file will be there.

## Build instructions:
1) Make a new directory build and cd into that directory
2) `cmake ..` (You may have to manually specifiy the location of waLBerla on your system)
3) Be sure to turn on code_gen and openmesh via e.g. ccmake, turn off tutorials/benchmarks/showcases
4) `make -j number_of_cores`
5) `git checkout .`  (To get the proper DynamicUBB.h header again)

## Run instructions
* You can run the application with a certain configuration file (`config.prm`) in the following ways:
1) `./aneurysm_model config.prm`
2) `mpirun -np $NUM_PROCESSES --bind-to core aneurysm_model config.prm` (with OpenMPI)

## Tips for running
1) Every heart beat cycle is 1 second physical time
2) The number of timesteps depends on the physical `dx`
3) `dt` scales as `O(dx^3)`, so the number of timesteps per heart beat cycle grows very quickly with finer meshes
4) Disable OpenMP for better performance (`export OMP_NUM_THREADS=1` also works)
5) The output files quickly become huge, the `resolution` variable can be used to control this

You can use Paraview to look at the results.

## Thanks
* A massive thank you to Markus Muhr (TUM) and Christoph Schwarzmeier (FAU) for their suggestions and improvements!
