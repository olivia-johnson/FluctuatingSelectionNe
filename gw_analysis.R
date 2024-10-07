## whole genome analysis

## Load packages and install any missing from the user's R library
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, tidyverse, ggpubr)

setwd("~/FluctuatingSelectionNe/")

# Vectors of initial loci number, epistasis parameter, target amplitude of fluctuation, and final loci number
loci=c(148, 127, 335,37,18)
epistasis=c(11, 15, 9.5, 2, 4.5)
amp_val=c(0.05,0.08, 0.02, 0.04, 0.35)
loci_val=c(111,90,264,27,9)

# create empty data object
af=NULL # allele frequency data set
ne_all=NULL # Ne data set

# for loop to compile simulation output
for (j in 1:length(loci)){ #loop through parameters to read in and compile files
  L=loci[j] # set initial loci number
  y=epistasis[j] # set y value
  amp=amp_val[j] # set amplitude
  Nl=loci_val[j] # set final loci number
  
  lab=paste0("l = ",L,", y = ", y) # set label of initial loci number and y value
  # read in replicates and add columns of the following constants
  for (rep in 1:10){ 
    # read in allele frequency file of no fitness simulation
    std_nofit=fread(paste0("timeseries_al_freq_ne_nofit_", rep,".txt")) 
    std_nofit[, sim_type:="Neutral"] # set simulation type
    std_nofit[, rep:=rep] # record replicate number
    std_nofit[, label:=lab] # record label
    std_nofit[, p_amp:=amp] # record expected amplitude
    std_nofit[, p_N:=Nl] # record final loci number
    
    # read in allele frequency file of no fitness capping simulation
    cap_nofit=fread(paste0("offcap_al_freq_ne_nofit_", rep,".txt"), fill=TRUE, skip="Gen") 
    cap_nofit[, sim_type:="Capped Neutral"] # set simulation type
    cap_nofit[, rep:=rep] # record replicate number
    cap_nofit[, label:=lab] # record label
    cap_nofit[, p_amp:=amp] # record expected amplitude
    cap_nofit[, p_N:=Nl] # record final loci number
    
    # read in Ne file of no fitness simulation
    n_ne=fread(paste0("timeseries_ne_nofit_", rep,".txt")) 
    n_ne[, sim_type:="Neutral"] # set simulation type
    n_ne[, rep:=rep] # record replicate number
    n_ne[, label:=lab] # record label
    n_ne[, p_amp:=amp] # record expected amplitude
    n_ne[, p_N:=Nl] # record final loci number
    
    # read in Ne file of no fitness capping simulation
    off_neu=fread(paste0("offcap_ne_nofit_", rep,".txt")) 
    off_neu[, sim_type:="Capped Neutral"] # set simulation type
    off_neu[, rep:=rep] # record replicate number
    off_neu[, label:=lab] # record label
    off_neu[, p_amp:=amp] # record expected amplitude
    off_neu[, p_N:=Nl] # record final loci number
    
    # read in allele frequency file of standard simulation
    std=fread(paste0("timeseries_al_freq_ne_", L,"_", y,"_", rep,".txt"))
    std[, sim_type:="Fluctuating"] # set simulation type
    std[, rep:=rep] # record replicate number
    std[, label:=lab] # record label
    std[, p_amp:=amp] # record expected amplitude
    std[, p_N:=Nl] # record final loci number
    
    # read in allele frequency file of capping simulation
    cap=fread(paste0("offcap_al_freq_ne_", L,"_", y,"_", rep,".txt"), fill=TRUE, skip="Gen")
    cap[, sim_type:="Capped Fluctuating"] # set simulation type
    cap[, rep:=rep] # record replicate number
    cap[, label:=lab] # record label
    cap[, p_amp:=amp] # record expected amplitude
    cap[, p_N:=Nl] # record final loci number
    
    # read in Ne file of standard simulation
    fs_ne=fread(paste0("timeseries_ne_",L,"_",y,"_", rep,".txt"))
    fs_ne[, sim_type:="Fluctuating"] # set simulation type
    fs_ne[, rep:=rep] # record replicate number
    fs_ne[, label:=lab] # record label
    fs_ne[, p_amp:=amp] # record expected amplitude
    fs_ne[, p_N:=Nl] # record final loci number
    
    # read in Ne file of capped simulation
    off_ne=fread(paste0("offcap_ne_",L,"_",y,"_", rep,".txt"))
    off_ne[, sim_type:="Capped Fluctuating"] # set simulation type
    off_ne[, rep:=rep] # record replicate number
    off_ne[, label:=lab] # record label
    off_ne[, p_amp:=amp] # record expected amplitude
    off_ne[, p_N:=Nl] # record final loci number
    
    af = rbind(af,std, cap, std_nofit, cap_nofit) # bind allele frequency data sets
    
    ne_all=rbind(ne_all, n_ne, off_ne, fs_ne, off_neu)} # bind Ne data sets
  
}


af[, year:=Gen%/%10, by=c("Gen")] # add value of year/seasonal cycle
af[, Freq.bin:="Segregating"] # set all alleles to segregating
af[mut_freq==0.000000, Freq.bin:="Fixed_Winter"] # relabel those fixed for winter allele
af[mut_freq==1.000000, Freq.bin:="Fixed_Summer"] # relabel those fixed for summer allele
af[Freq.bin=="Segregating",n_seg:=.N, by = c("Gen", "sim_type", "rep", "label")] # count the number of segregating loci
af[Freq.bin=="Segregating",fc_year:=max(mut_freq)-min(mut_freq), by=c("mut_id", "year", "rep", "sim_type", "label")] #calculate fequency change over the year of segregating alleles

# calculate mean, max and median allele frequency change
mean_amp=af[Gen>=10000 & Freq.bin=="Segregating", .(mean=mean(fc_year), max=max(fc_year), median=median(fc_year)), by=c("sim_type", "label", "rep")]
mean_amp[, amp:=round(mean*100, digits = 1)] # make amplitude a percentage
# identify loci segregating at the final time point
final_loci=af[Gen==10010 & Freq.bin=="Segregating", unique(n_seg), by=c("rep", "sim_type", "label")]
# calculate the mean number if segregating loci  for each parameter combination and simulation type
final_loci[, .(mean(V1)), by=c("label", "sim_type")]
# calculate the difference between expected and actual number of final segregating loci
af[, n_diff:=n_seg-p_N]
# calculate the mean frequency change for each replicate
af[Freq.bin=="Segregating", rep_amp:=mean(fc_year), by=c("rep", "sim_type", "label", "year")]
# compare the mean frequency fluctuation with the expected allele frequency fluctuation
af[, amp_diff:=rep_amp-p_amp]
# summarise the difference between expected and actual values for the simulations (i.e. difference in no. seg alleles, amplitude etc)
af[Gen==10009 & sim_type=="Fluctuating"&Freq.bin=="Segregating", .(mean_diff=mean(n_diff), max_diff=max(n_diff), mean_amp_diff=mean(abs(amp_diff)), max_amp_diff=max(abs(amp_diff))), by="label"]
# overall difference between the derived and actual parameters
af[Gen==10009 & sim_type=="Fluctuating"&Freq.bin=="Segregating", .(mean_diff=mean(n_diff), max_diff=max(n_diff), mean_amp_diff=mean(abs(amp_diff)), max_amp_diff=max(abs(amp_diff)))]

# Extract Ne values and relevant metadata
ne_data=ne_all[,.(Time, linked_Ne, unlinked_Ne, sim_type, rep, label)]
# transform data set so Ne all in one column
ne_data=melt(ne_data, id.vars = c("Time", "sim_type", "rep", "label"), value="ne", variable.name = "linkage")
## reformat linkage label
ne_data[,linkage:= tstrsplit(linkage, "_", keep=1)]
# calculate mean Ne for each parameter set and simulation type
ne_data[, mean_ne:=mean(ne), by=c("Time", "linkage", "sim_type", "label")]
# calculate the mean reduction in Ne (500,000 individuals in simulated populations)
ne_data[, reduction:=mean_ne/500000-1, by=c("Time", "linkage", "sim_type", "label")]
# calculate the mean reduction in Ne for each replicate(500,000 individuals in simulated populations)
ne_data[, rep_reduction:=ne/500000-1, by=c("Time", "linkage", "sim_type", "label", "rep")]

# Figure 4 - Relative reduction in effective population size due to fluctuating selection

# create amplitude annotations
relne.labs=mean_amp[sim_type=="Fluctuating", .( amp=round(mean(amp),1)), by="label"]

rel_ne=ggplot(ne_data[linkage=="unlinked" & sim_type!="Capped Fluctuating" & sim_type!="Capped Neutral"], aes(x=Time))+
  geom_line(aes(y=reduction, group=sim_type, col=sim_type), linewidth=1)+
  geom_hline(yintercept = 0, linetype="dashed")+
  theme_bw()+coord_cartesian(y=c(-1, .1))+
  labs(x="Generations", y="Relative Reduction in Ne", col="Selection type") + 
  facet_wrap(~label, nrow=1)+
  geom_text(data= relne.labs,mapping = aes(x=6000, y=-0.95, label = paste("Mean Amp =",amp, "%")), size=3)
ggsave(paste0("plots/Figure4.pdf"), rel_ne, width=10, height=3)
ggsave(paste0("plots/Figure4.jpg"), rel_ne, width=10, height=3)
  
ne_data[, id:=paste0(sim_type, "", rep)] # add id to separate replicate lines
rel_ne_var= ggplot(ne_data[linkage=="unlinked" & sim_type!="Capped Fluctuating" & sim_type!="Capped Neutral"  ], aes(x=Time))+
  geom_line(aes(y=rep_reduction, group=id, col=sim_type), linewidth=0.3, alpha=0.5)+
  geom_line(aes(y=reduction, group=id, col=sim_type), linewidth=1)+
  geom_hline(yintercept = 0, linetype="dashed")+
  theme_bw()+coord_cartesian(y=c(-1, .1))+
  labs(x="Generations", y="Relative Reduction in Ne", col="Selection type") + 
  facet_wrap(~label, nrow=1)+
  geom_text(data= relne.labs,mapping = aes(x=6000, y=-0.95, label = paste("Mean Amp =",amp, "%")), size=3)
ggsave(paste0("plots/FigureS2.pdf"), rel_ne, width=10, height=3)
ggsave(paste0("plots/FigureS2.jpg"), rel_ne, width=10, height=3)

# Figure 5 - Effective population size across the seasonal cycle.

 # compile seasonal ouput
# create empty data objects
s_af=NULL # seasonal allele frequencies
sne_all=NULL # seasonal Ne

for (j in 1:length(loci)){ # loop through parameter combinations
  L=loci[j] # set initial loci number
  y=epistasis[j] # set epistasis value
  amp=amp_val[j] # set amplitude of fluctuation
  Nl=loci_val[j] # set final loci number
  
  lab=paste0("l = ",L,", y = ", y) # create label
  for (rep in 1:10){ # for each replicate
    # read in seasonal allele frequency data
    seas=fread(paste0("season_al_freq_ne_", L,"_", y,"_", rep,".txt"))
    seas[, sim_type:="Fluctuating"] # set simulation type
    seas[, rep:=rep] # set replicate number
    seas[, label:=lab] # set label
    seas[, p_amp:=amp] # set amplitude
    seas[, p_N:=Nl] # set final loci number
    
    # read in seasonal Ne data
    seas_ne=fread(paste0("season_ne_",L,"_",y,"_", rep,".txt"))
    seas_ne[, sim_type:="Fluctuating"] # set simulation type
    seas_ne[, rep:=rep] # set replicate number
    seas_ne[, label:=lab] # set label
    seas_ne[, p_amp:=amp] # set amplitude
    seas_ne[, p_N:=Nl] # set final loci number
    
    s_af = rbind(s_af,seas) # bind to seasonal allele frequency data set
    
    sne_all=rbind(sne_all,seas_ne)} # bind to seasonal Ne data set
}

s_af[, year:=Gen%/%10, by=c("Gen")] # define year/seasonal cycle
s_af[, Freq.bin:="Segregating"] # set all alleles to segregating
s_af[mut_freq==0.000000, Freq.bin:="Fixed_Winter"] # relabel alleles fixed for winter
s_af[mut_freq==1.000000, Freq.bin:="Fixed_Summer"] # relabel alleles fixed for summer
s_af[Freq.bin=="Segregating",n_seg:=.N, by = c("Gen", "sim_type", "rep", "label")] # count number of segregating alleles
# calculate frequency change over the year
s_af[Freq.bin=="Segregating",fc_year:=max(mut_freq)-min(mut_freq), by=c("mut_id", "year", "rep", "sim_type", "label")] 
# calculate the mean amplitude for each replicate
s_af[Freq.bin=="Segregating", rep_amp:=mean(fc_year), by=c("rep", "sim_type", "label", "year")]
s_af[, amp_diff:=rep_amp-p_amp] # calculate difference between replciate amplitude and expected amplitude
# isolate data required for plots
sne_data=sne_all[,.(Time, linked_Ne, unlinked_Ne, sim_type, rep, label)]
# transform data table so Ne is in one column
sne_data=melt(sne_data, id.vars = c("Time", "sim_type", "rep", "label"), value="ne", variable.name = "linkage")
# fix linkage labels
sne_data[,linkage:= tstrsplit(linkage, "_", keep=1)]
# calculate mean Ne for parameter set and simulation type
sne_data[, mean_ne:=mean(ne), by=c("Time", "linkage", "sim_type", "label")]
# calculate mean overal reduction in Ne
sne_data[, reduction:=mean_ne/500000-1, by=c("Time", "linkage", "sim_type", "label")]
# calculate reduction in Ne for each replicate
sne_data[, rep_reduction:=ne/500000-1, by=c("Time", "linkage", "sim_type", "label", "rep")]

sne_data[, id:=paste0(label,rep), by=c("label", "rep")] # create label to separate replicate lines in plot
s_ne=ggplot()+geom_vline(xintercept = 2505, col="blue", linetype="dotted")+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_line(data=sne_data[linkage=="unlinked"],aes(x=Time,y=reduction, group=label), linewidth=1)+
  geom_line(data=sne_data[linkage=="unlinked"],aes(x=Time,y=rep_reduction, group=id, col=label), linewidth=0.5, alpha=0.5)+
  theme_bw()+coord_cartesian(x=c(2500,2510))+ scale_x_continuous(breaks=c(2500,2505,2510))+ 
  labs(x="Generations", y="Relative Reduction in Ne", col="Simulation Parameters") + facet_wrap(~label,nrow=1)+theme(legend.position = 'none')
ggsave(paste0("plots/FigureS3.pdf"), s_ne, width=10, height=3)
ggsave(paste0("plots/FigureS3.jpg"), s_ne, width=10, height=3)

s_ne_mean=ggplot()+geom_vline(xintercept = 2505, col="blue", linetype="dotted")+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_line(data=sne_data[linkage=="unlinked"],aes(x=Time,y=reduction, group=label, col=label), linewidth=1)+
  theme_bw()+coord_cartesian(x=c(2500,2510))+ scale_x_continuous(breaks=c(2500,2502,2504,2506,2508,2510))+ 
  labs(x="Generations", y="Relative Reduction in Ne", col="Selection Parameters") 

ggsave(paste0("plots/Figure5.pdf"), s_ne_mean, width=5, height=2.5)
ggsave(paste0("plots/Figure5.jpg"), s_ne_mean, width=5, height=2.5)


# Figure 6 - Relationship between reduction in effective population size and aspects of seasonal fluctuation (amplitude).
# merge Ne and amplitude data of standard simulations
ne_amp=merge(ne_data[Time==10010 & linkage=="unlinked"], mean_amp, by=c("sim_type", "label", "rep"))
# plot mean, median and maximum amplitude vs amplitude
mean=ggplot(ne_amp[sim_type=="Fluctuating" ])+
  geom_point(aes(x=mean, y=rep_reduction, col=label)) + labs(x="Mean Amplitude", y="Reduction in Ne", col="Parameters") + theme_bw()
med=ggplot(ne_amp[sim_type=="Fluctuating"])+
  geom_point(aes(x=median, y=rep_reduction, col=label))+ labs(x="Median Amplitude", col="Parameters")+ theme_bw()+theme(axis.title.y = element_blank())
max=ggplot(ne_amp[sim_type=="Fluctuating"])+
  geom_point(aes(x=max, y=rep_reduction, col=label)) + 
  labs(x="Maximum Amplitude", col="Parameters")+ theme_bw() +theme(axis.title.y = element_blank())

ne_amp_plot=ggarrange(mean, med, max,common.legend = TRUE, legend = "right", nrow=1 )
ggsave(paste0("plots/Figure6.pdf"), ne_amp_plot, width=8, height=3)
ggsave(paste0("plots/Figure6.jpg"), ne_amp_plot, width=8, height=3)

# Linear regression to predict reduction in Ne from empirical values (calculated in empiricalMaxAmp.R)
max.mod=lm(rep_reduction~max, ne_amp[sim_type=="Fluctuating"])
summary(max.mod)

machado_max=data.table(max=0.45) # value calculated in in empiricalMaxAmp.R
mach_ne=max.mod %>% predict(machado_max) # predict reduction in Ne in empirical populations from Machado et al. 2021
bergland_max=data.table(max=0.37) # value calculated in in empiricalMaxAmp.R
berg_ne=max.mod %>% predict(bergland_max) # predict reduction in Ne in empirical populations from Bergland et al. 2014

# Figure 7 - Mean effect size and mean dominance coefficient of seasonal loci segregating at the end of simulations in relation to their amplitude of fluctuation.

effectsize=NULL # create empty data object to store effect size and dominance information

for (j in 1:length(loci)){ # loop through parameter sets
  L=loci[j] # set initial loci number
  y=epistasis[j] # set y value
    for (rep in 1:10){ # loop through replicates
    muts=fread(paste0("seasonal_mutations_",L,"_",y,"_", rep,".txt")) # read in mutation information
    muts[, rep:=rep] # set replicate number
    muts[, label:=paste0("l = ",L,", y = ", y)] # define label
    muts[, mut_id:=0:(L-1)] # give each mutation, individual id
    effectsize=rbind(effectsize, muts)}} # bind to data set
# isolate allele frequency information for fluctuating and capped fluctuating sims
final_gen=af[Gen==10009  & sim_type!="Neutral" & sim_type!="Capped Neutral"] 
# merge allele frequencies and mutation information
final_gen=left_join(final_gen, effectsize, by=c("rep", "label", "mut_id"), multiple="all")
final_gen[, mean_fx:=(s_fx+w_fx)/2] # calculate mean effect size across seasons
final_gen[, mean_d:=(s_d+w_d)/2] # calculate mean dominance across seasons

fx_amp=ggplot(final_gen[sim_type=="Capped Fluctuating"], aes(x=mean_fx, y=fc_year, col=mean_d))+
  geom_point(alpha=0.8) + facet_wrap(~label, nrow=1)+scale_x_log10()+scale_y_log10() + 
  scale_color_gradientn(trans="log10", na.value = "white",colours=c("orange", "darkmagenta", "navyblue"))+
  labs(x = "Mean Effect size", y="Amplitude", col= "Mean Dominance")+theme_bw()

d_amp=ggplot(final_gen[ sim_type=="Capped Fluctuating"], aes(x=mean_d, y=fc_year, col=mean_fx))+
  geom_point(alpha=0.8) + facet_wrap(~label, nrow=1)+scale_x_log10()+scale_y_log10() + 
  scale_color_gradientn(trans="log10", na.value = "white",colours=c("orange", "darkmagenta", "navyblue"))+
  labs(x = "Mean Dominance Coefficient", y="Amplitude", col="Mean Effect size") +theme_bw()

f_d=ggarrange(fx_amp, d_amp, labels = c("A", "B"), ncol=1)#,common.legend = TRUE, legend = "right")
ggsave(paste0("plots/Figure7.pdf"), f_d, width=10, height=6)
ggsave(paste0("plots/Figure7.jpg"), f_d, width=10, height=6)


# Figure 8 - Reduction in effective population size for varying census size (Nc).
# create empty data objects
nc_af=NULL # allele frequency data object
nc_all=NULL # ne data object
 # compile data
    for (rep in 1:10){ # for each replicate (of 10)
    # read in allele frequencies for standard simulation with Nc 500,000
    full=fread(paste0("timeseries_al_freq_ne_127_15_",rep,".txt"))
    full[, sim_type:="500000"] # define Nc
    full[, rep:=rep] # set replicate number
    full[, p_amp:=0.08] # set expected amplitude
    full[, p_N:=90] # set  expected final loci number
    
    # read in allele frequencies for simulation with Nc 10,000
    tt=fread(paste0("timeseries_al_freq_ne_127_15_10000_",rep,".txt"), fill=TRUE, skip="Gen")
    tt[, sim_type:="10000"] # define Nc
    tt[, rep:=rep] # set replicate number
    tt[, p_amp:=0.08] # set expected amplitude
    tt[, p_N:=90] # set  expected final loci number
    
    # read in allele frequencies for simulation with Nc 50,000
    ft=fread(paste0("timeseries_al_freq_ne_127_15_50000_",rep,".txt"))
    ft[, sim_type:="50000"] # define Nc
    ft[, rep:=rep] # set replicate number
    ft[, p_amp:=0.08] # set expected amplitude
    ft[, p_N:=90] # set  expected final loci number
    
    # read in allele frequencies for simulation with Nc 200,000
    oht=fread(paste0("timeseries_al_freq_ne_127_15_100000_",rep,".txt"))
    oht[, sim_type:="100000"] # define Nc
    oht[, rep:=rep] # set replicate number
    oht[, p_amp:=0.08] # set expected amplitude
    oht[, p_N:=90] # set  expected final loci number
    
    # read in Ne for standard simulation with Nc 500,000
    full_ne=fread(paste0("timeseries_ne_127_15_",rep,".txt"))
    full_ne[, sim_type:="500000"] # define Nc
    full_ne[, rep:=rep] # set replicate number
    full_ne[, p_amp:=0.08] # set expected amplitude
    full_ne[, p_N:=90] # set  expected final loci number
    
    # read in Ne for  simulation with Nc 10,000
    tt_ne=fread(paste0("timeseries_ne_127_15_10000_",rep,".txt"))
    tt_ne[, sim_type:="10000"] # define Nc
    tt_ne[, rep:=rep] # set replicate number
    tt_ne[, p_amp:=0.08] # set expected amplitude
    tt_ne[, p_N:=90] # set  expected final loci number
    
    # read in Ne for  simulation with Nc 50,000
    ft_ne=fread(paste0("timeseries_ne_127_15_50000_",rep,".txt"))
    ft_ne[, sim_type:="50000"] # define Nc
    ft_ne[, rep:=rep] # set replicate number
    ft_ne[, p_amp:=0.08] # set expected amplitude
    ft_ne[, p_N:=90]  # set  expected final loci number
    
    # read in Ne for  simulation with Nc 100,000
    oht_ne=fread(paste0("timeseries_ne_127_15_100000_",rep,".txt"))
    oht_ne[, sim_type:="100000"] # define Nc
    oht_ne[, rep:=rep] # set replicate number
    oht_ne[, p_amp:=0.08] # set expected amplitude
    oht_ne[, p_N:=90] # set  expected final loci number
    
    nc_af = rbind(nc_af,full, tt, ft, oht) # bind allele frequency data sets
    
    nc_all=rbind(full_ne, tt_ne, ft_ne, oht_ne, nc_all)} # bind Ne datasets
  
nc_af[, year:=Gen%/%10, by=c("Gen")] # define year/seasonal cycle
nc_af[, Freq.bin:="Segregating"] # set all alleles to segregating
nc_af[mut_freq==0.000000, Freq.bin:="Fixed_Winter"] # relabel allele fixed for winter
nc_af[mut_freq==1.000000, Freq.bin:="Fixed_Summer"] # relabel alleles fixed for summer
nc_af[Freq.bin=="Segregating",n_seg:=.N, by = c("Gen", "sim_type", "rep")] # count number of segregating alleles
nc_af[Freq.bin=="Segregating",fc_year:=max(mut_freq)-min(mut_freq), by=c("mut_id", "year", "rep", "sim_type")] # calculate frequency change
# calculate mean, maximum and median fluctuation
mean_ncamp=nc_af[Gen>=10000 & Freq.bin=="Segregating", .(mean=mean(fc_year), max=max(fc_year), median=median(fc_year)), by=c("sim_type", "rep")]
mean_ncamp[, amp:=round(mean*100, digits = 1)] # reformat amplitude as percentage
nc_af[, n_diff:=n_seg-p_N] # calculate difference in expected and actual number of segregating alleles
# calculate mean fluctuation for each replicate
nc_af[Freq.bin=="Segregating", rep_amp:=mean(fc_year), by=c("rep", "sim_type",  "year")]
# calculate difference in expected and actual frequency of fluctuation
nc_af[, amp_diff:=rep_amp-p_amp]
# subset Ne data to what is reuqire for plotting
nc_data=nc_all[,.(Time, linked_Ne, unlinked_Ne, sim_type, rep)]
# transform dataset so Ne is in one column
nc_data=melt(nc_data, id.vars = c("Time", "sim_type", "rep"), value="ne", variable.name = "linkage")
#reformat linkage labels
nc_data[,linkage:= tstrsplit(linkage, "_", keep=1)]
# calculate mean Ne per Nc
nc_data[, mean_ne:=mean(ne), by=c("Time", "linkage", "sim_type")]
# calculate the mean reduction in Ne for each Nc
nc_data[, reduction:=mean_ne/as.numeric(sim_type)-1, by=c("Time", "linkage", "sim_type")]
# calculate reduction in Ne for each replicate
nc_data[, rep_reduction:=ne/as.numeric(sim_type)-1, by=c("Time", "linkage", "sim_type",  "rep")]
# create amplitude labels for plots
relnc.labs=mean_ncamp[, .(amp=mean(mean)), by="sim_type"]
# create labels for facet strips
nc_data[, label:=paste("Nc =", sim_type), by="sim_type"]
# create labels for facet strips
relnc.labs[, label:=paste("Nc =", sim_type), by="sim_type"]


rel_nc=ggplot(nc_data[linkage=="unlinked"], aes(x=Time))+
  geom_line(aes(y=rep_reduction, group=rep,col=sim_type), linewidth=0.3, alpha=0.5)+
  geom_line(aes(y=reduction, col=sim_type), linewidth=1)+
  geom_hline(yintercept = 0, linetype="dashed")+
  theme_bw()+coord_cartesian(y=c(-1, .1))+
  labs(x="Generations", y="Relative Reduction in Ne", col="Selection type") + 
  facet_wrap(~label, nrow=1)+ theme(legend.position="none")+
  geom_text(data= relnc.labs,mapping = aes(x=6000, y=-0.95, label = paste("Mean Amp =",signif(amp,2), "%")), size=3)
ggsave(paste0("plots/Figure8.pdf"), rel_nc, width=10, height=3)
ggsave(paste0("plots/Figure 8.jpg"), rel_nc, width=10, height=3)


# Figure 9 - Reduction in effective population size with and without offspring capping.
# create labels for plots of mean amplitude for capped fluctuating simulations
caprelne.labs=mean_amp[sim_type=="Capped Fluctuating", .(amp=round(mean(amp),1)), by="label"]

rel_ne=ggplot(ne_data[linkage=="unlinked"], aes(x=Time))+
  geom_line(aes(y=reduction, group=id, col=sim_type), linewidth=1)+
  geom_hline(yintercept = 0, linetype="dashed")+
  theme_bw()+coord_cartesian(y=c(-1, .1))+
  labs(x="Generations", y="Relative Reduction in Ne", col="Selection type") + 
  facet_wrap(~label, nrow=1)+
  geom_text(data= caprelne.labs,mapping = aes(x=6000, y=-0.95, label = paste("Mean Amp =",amp, "%")), size=3)
ggsave(paste0("plots/Figure9.pdf"), rel_ne, width=10, height=3)
ggsave(paste0("plots/Figure9.jpg"), rel_ne, width=10, height=3)

# summarise capped vs standard amplitudes
capmean=mean_amp[sim_type=="Capped Fluctuating"| sim_type=="Fluctuating", .(mean=mean(mean), max=max(max), median=mean(median),ID=paste0(sim_type,label)), by=c("sim_type", "label")]

# Figure 10 - Reduction in effective population size due to fluctuating selection in a population with boom-bust demography.
# create empty data objects
fp.af=NULL # fluctuating population size allele frequency
fp.ne_all=NULL # fluctuating population size Ne

for (j in 1:length(loci)){ # loop through parameter sets
  L=loci[j] # set initial loci number
  y=epistasis[j] # set y value
  
  lab=paste0("FP_",L,"_", y) # define label in file names
  rlab=paste0("l = ",L,", y = ", y) # define labels for plots
  
  for (rep in 1:10){ # for each replicate
    # read in allele frequency file no fitness
    std_nofit=fread(paste0("timeseries_al_freq_ne_FP_nofit_", rep,".txt"))
    std_nofit[, sim_type:="Fluctuating Pop\nNeutral"] # define simulation type
    std_nofit[, rep:=rep] # define replicate number
    std_nofit[, label:=rlab] # add label for plotting
    
    # read in Ne file no fitness
    n_ne=fread(paste0("timeseries_ne_FP_nofit_", rep,".txt"))
    n_ne[, sim_type:="Fluctuating Pop\nNeutral"] # define simulation type
    n_ne[, rep:=rep] # define replicate number
    n_ne[, label:=rlab] # add label for plotting
    
    # read in allele frequency file fluctuating
    std=fread(paste0("timeseries_al_freq_ne_", lab,"_", rep,".txt"))
    std[, sim_type:="Fluctuating Pop\nFluctuating"] # define simulation type
    std[, rep:=rep] # define replicate number
    std[, label:=rlab] # add label for plotting
    
    # read in Ne file fluctuating
    fs_ne=fread(paste0("timeseries_ne_", lab,"_", rep,".txt"))
    fs_ne[, sim_type:="Fluctuating Pop\nFluctuating"] # define simulation type
    fs_ne[, rep:=rep] # define replicate number
    fs_ne[, label:=rlab] # add label for plotting
    
    fp.af = rbind(fp.af,std, std_nofit) # bind alelle frequency objectes
    
    fp.ne_all=rbind(fp.ne_all, n_ne,  fs_ne)} # bind ne objects
}

fp.af[, year:=Gen%/%10, by=c("Gen")] # define year / seasonal cycle
fp.af[, Freq.bin:="Segregating"] # set all alleles to segregating
fp.af[mut_freq==0.000000, Freq.bin:="Fixed_Winter"] # relabel alleles fixed for winter 
fp.af[mut_freq==1.000000, Freq.bin:="Fixed_Summer"] # relabel allelles fixed for summer
fp.af[Freq.bin=="Segregating",n_seg:=.N, by = c("Gen", "sim_type", "rep", "label")] # count number of segregating loci
# calculate frequency fluctuation
fp.af[Freq.bin=="Segregating",fc_year:=max(mut_freq)-min(mut_freq), by=c("mut_id", "year", "rep", "sim_type", "label")]
# calculate mean fluctuating for segregating alleles
fp_mean_amp=fp.af[Gen>=10000 & Freq.bin=="Segregating", mean(fc_year), by=c("sim_type", "label")]
fp_mean_amp[, amp:=V1*100] # convert to percentage
h.mean=10/((1/5e5)*5+(1/5e4)*5) # calculate harmonic mean population size
# subset datasheet to that required for plotting
fp.ne_data=fp.ne_all[,.(Time, linked_Ne, unlinked_Ne, sim_type, rep, label)]
# transform data so Ne in one column
fp.ne_data=melt(fp.ne_data, id.vars = c("Time", "sim_type", "rep", "label"), value="ne", variable.name = "linkage")
# reformat linkage label
fp.ne_data[,linkage:= tstrsplit(linkage, "_", keep=1)]
# calculate mean Ne
fp.ne_data[, mean_ne:=mean(ne), by=c("Time", "linkage", "sim_type", "label")]
# calculate the reduction in Ne relative to the harminc mean population size
fp.ne_data[, reduction:=mean(ne)/h.mean-1, by=c("Time", "linkage", "sim_type", "label")]
# create labels of mean amplitude for fluctuating population size
fprelne.labs=fp_mean_amp[sim_type=="Fluctuating Pop\nFluctuating", .(label, amp)]
cp=ne_data[linkage=="unlinked" & sim_type=="Fluctuating" | sim_type=="Neutral"] # extract constant population size data
cp[,(c( "rep_reduction", "id")):=NULL] # remove columns to be able to  ind for fluctuating population size data
# fix labels names
cp[, sim_type:=ifelse(sim_type=="Fluctuating", "Constant Pop\nFluctuating", "Constant Pop\nNeutral"), by="sim_type"]
#bind constant and fluctuation population size data for plot
fp.ne_data=rbind(fp.ne_data, cp)

rel_ne=ggplot()+
  geom_line(data=fp.ne_data[linkage=="unlinked"],aes(x=Time,y=reduction, group=sim_type, col=sim_type), linewidth=1)+
  geom_hline(yintercept = 0, linetype="dashed")+
  theme_bw()+
  labs(x="Generations", y="Relative Reduction in Ne", col="Selection type") +
  facet_wrap(~label, nrow=1)+ 
  geom_text(data= fprelne.labs,mapping = aes(x=6000, y=-0.95, label = paste("Mean Amp =",round(amp, digits=2), "%")), size=3)
ggsave(paste0("plots/Figure10.pdf"), rel_ne, width=10, height=3)
ggsave(paste0("plots/Figure10.jpg"), rel_ne, width=10, height=3)

# calculate Ne using harminc population size and reduction in Ne for different parameter sets
h.mean*(1+unique(fp.ne_data[Time==10010 & sim_type=="Fluctuating Pop\nFluctuating" & linkage=="unlinked" &label=="l = 127, y = 15", reduction]))
h.mean*(1+unique(fp.ne_data[Time==10010 & sim_type=="Fluctuating Pop\nFluctuating" & linkage=="unlinked" &label=="l = 148, y = 11", reduction]))
h.mean*(1+unique(fp.ne_data[Time==10010 & sim_type=="Fluctuating Pop\nFluctuating" & linkage=="unlinked" &label=="l = 37, y = 2", reduction]))
h.mean*(1+unique(fp.ne_data[Time==10010 & sim_type=="Fluctuating Pop\nFluctuating" & linkage=="unlinked" &label=="l = 18, y = 4.5", reduction]))
h.mean*(1+unique(fp.ne_data[Time==10010 & sim_type=="Fluctuating Pop\nFluctuating" & linkage=="unlinked" &label=="l = 335, y = 9.5", reduction]))
