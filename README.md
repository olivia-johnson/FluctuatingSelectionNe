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
slim -d results_dir=${results_dir} -d fit=${fitness} -d  L=${oci} -d y=${epistasis} -d rep=${rep} ~/FluctuatingSelectionNe/timeseries_Ne_short.slim
```
