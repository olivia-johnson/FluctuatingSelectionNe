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

### Genome-wide simulations

To run simulations use the following code 

```ruby
slim -d results_dir=${results_dir} -d fit=${fitness} -d  L=${loci} -d y=${epistasis} -d rep=${rep} timeseries_Ne_short.slim
```
Where results_dir is the path where you want output to be written to, fitness is either 0 (to turn off fitness function) or 1 (to engage fluctuating fitness model), loci is the number of initial seasonal loci you want drawn onto the genome, epistasis is the parameter _y_, and rep is the replicate number.
