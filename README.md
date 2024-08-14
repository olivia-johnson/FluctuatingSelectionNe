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

Once the parameter sets have been define, the simulation can be run using the following code,

```ruby
python chp4_allele.py ${params} ${replicate number} ${results_dir}
```
This will run [chp4_allele.py](chp4_allele.py), which will read in and parse parameters to the [SLiM simulation](witt_complex_allele.slim).   

### Genome-wide simulations

To run simulations use the following code 

```ruby
slim -d results_dir=${results_dir} -d fit=${fitness} -d  L=${loci} -d y=${epistasis} -d rep=${rep} timeseries_Ne_short.slim
```
Where results_dir is the path where you want output to be written to, fitness is either 0 (to turn off fitness function) or 1 (to engage fluctuating fitness model), loci is the number of initial seasonal loci you want drawn onto the genome, epistasis is the parameter _y_, and rep is the replicate number.

> [!IMPORTANT]
> You must first run general genome-wide simulations with fitness = 1 to determine the seasonal dominance and effect size of each locus. These values are used across simulation types (i.e. capping) but are only generated in [timeseries_Ne_short.slim](timeseries_Ne_short.slim) simulations with fitness acting. 
