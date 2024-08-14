path='~/FluctuatingSelectionNe/'
import os
import sys
import msprime
import pyslim
import numpy as np
import yaml
import tskit
import pandas as pd
import time

params=sys.argv[1] ## define parameter set 
sim_run = sys.argv[2] ## define simulation replicate number
tmpdir =str(sys.argv[3]) ## define directory to write output files

def simulate_alleles(tmpdir, group, sim_run, s_pop, w_pop, l, y, rGen, fitness_on, sum_gen, win_gen, path):
        genomeSize = l ## genome size is the length of the number of loci

        ## FORWARD SIMULATION
        start_time = time.time()
        tmpdir_call = "tmpdir='" + str(tmpdir)+ "'" 
        ## command to run SLiM simulation
        cmd = 'slim -d "' +str(tmpdir_call)+ '" -d fit='+ str(fitness_on)+" -d group=" + str(group) + " -d g_s=" + str(sum_gen)+" -d g_w=" + str(win_gen)+" -d sim_run=" + str(sim_run) + " -d GenomeSize=" + str(int(genomeSize)) + " -d L=" + str(l)+ " -d n_s=" + str(int(s_pop)) + " -d n_w=" + str(int(w_pop)) + " -d y=" + str(y) + " -d rGen="+ str(rGen) +path+"/witt_complex_allele.slim"
        print(cmd)
        os.system(cmd) ## runn command
        print("Simulations took ", (time.time()-start_time) ,  " seconds")


####  READ IN PARAMETERS
    # load in parameter file
with open('{0}parameters/{1}.txt'.format(path,params), 'r') as f: ## load in parameter file
    parameters = yaml.load(f, Loader=yaml.FullLoader)

    #set parameters from file
s_pop = int(parameters["s_pop"])
w_pop = int(parameters["w_pop"])
l = int(parameters["l"])
y = parameters["y"]
rGen=int(parameters["rGen"])
fitness_on = parameters["fitness_on"]
sum_gen = int(parameters["sum_gen"])
win_gen = int(parameters["win_gen"])
group=parameters["group"]

start_time = time.time()

####  SIMULATE  ALLELE FREQUENCIES - run function define above
simulate_alleles(tmpdir, group, sim_run, s_pop, w_pop, l, y, rGen, fitness_on, sum_gen, win_gen, path)

print("Time = ", (time.time() - start_time))