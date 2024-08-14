# The effect of fluctuating selection on effective population size (N<sub>e</sub>)
This repository contains the code to replicate the analysis conducted in the abovenamed manuscript. 

## Directions

### Setting up the environment
The [yaml](ne_env.yml) file can be used to create a conda environment with the programs required to conduct the simulations and analysis performed in this study.

Once this repository has been cloned, the environment can be created using the following code

```ruby
conda env create -f ne_env.yml

conda activate ne_env
```

### Conducting simulations of fluctuating allele frequencies
Parameters for allele frequency simulations are recorded in a text file labelled [group_x.txt](group_x.txt), where x is replaced by a unique numeric identifier for each parameter set. This file can be used as a template for the parameter file. The values we used for variables are included in the comments of this file in brackets, with constant values already included in the parameter file. Remove the text in quotation marks before running the simulations.

Once the parameter sets have been defined, the simulation can be run using the following code,

```ruby
python runAllele.py ${group_x} ${rep} ${results_dir}
```
group_x is the parameter set label, rep is the replicate number and results_dir is the path where you want output to be written to.

This will run [runAllele.py](runAllele.py), which will read in and parse parameters to the [SLiM simulation](wittmann_allele.slim). 

**You may need to change the path at the start of [runAllele.py](runAllele.py)**

Each simulation will produce a file named al_freq_group_x_rep.txt

> [!NOTE]
> These simulations (particularly with higher loci numbers) are time and memory-intensive.

### Genome-wide simulations

To run constant population size simulations use the following code, 

```ruby
slim -d results_dir=${results_dir} -d fit=${fitness} -d  L=${loci} -d y=${epistasis} -d rep=${rep} gw_sim.slim
```
results_dir is the path where you want output to be written to, fitness is either 0 (to turn off fitness function) or 1 (to engage fluctuating fitness model), loci is the number of initial seasonal loci you want drawn onto the genome, epistasis is the parameter _y_, and rep is the replicate number.

For fluctuating population size, you also need to specify the winter population size (n_w),

``` ruby
slim -d results_dir=${results_dir} -d n_w=50000 -d fit=${fitness} -d  L=${loci} -d y=${epistasis} -d rep=${rep} gw_sim.slim
```
To run simulations with offspring capping or sample only across a single seasonal cycle, replace ``` gw_sim.slim ``` with either ``` capping.slim ``` or ``` seasonal_ne.slim ``` in the above code.

> [!IMPORTANT]
> You must first run constant a population size genome-wide simulation with fitness = 1 to determine the seasonal dominance and effect size of each locus. These values are used across simulation types (i.e. [capping](capping.slim) and [seasonal](seasonal_ne.slim)) but are only generated in [gw_sim.slim](gw_sim.slim) simulations with fitness acting.
>
> Each replicate of each parameter set (loci number and _y_ value) will have a different set of seasonal mutations that will be used across simulation types (i.e. seasonal mutations are constant for replicate 1 across the different simulations).
