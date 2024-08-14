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

To run simulations use the following code 

```ruby
slim -d results_dir=${results_directory} -d fit=${fitness (0 for no fitness or 1 to enact fluctuating fitness model)} -d  L=${number of initial loci} -d y=${epistasis parameter} -d rep=${replicate number} ~/FluctuatingSelectionNe/timeseries_Ne_short.slim
```
