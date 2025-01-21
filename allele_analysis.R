 ### Examine allele frequency fluctuations and generate multiple regression models###

## Load packages and install any missing from the user's R library
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, ggplot2, plotly)

path="~/FluctuatingSelectionNe/" # set path 
setwd(path) # set working directory
groups=c(1:x) # parameter set IDs

## Compile allele frequency files into one data set with parameter values
alfreq_data = NULL 
for (g in groups){
  print(g)
  # list ferquency files
  f_list  <- list.files(path =paste0(path,"/group_", g, "/"),pattern ="al_freq_")
  # list parameter files
  parameters <- fread(file=paste0(path,"parameters/group_",g, ".txt"), sep = ":")
  setkey(parameters, V1)
  # read in parameters
  fiton=parameters["fitness_on", V2]
  if (fiton==0){
    epi="Drift"
  } else{
    epi = parameters["y", V2]
  }
  sum_gen =parameters["sum_gen", V2]
  win_gen =parameters["win_gen", V2] 
  if (sum_gen==win_gen){
    gen_s = "EG"
  } else{
    gen_s = "UG"
  }
  sum_pop =parameters["s_pop", V2]
  win_pop =parameters["w_pop", V2]
  if (sum_pop==win_pop){
    pop_s = "CP"
  } else{
    pop_s = "FP"
  }
  loci =parameters["l", V2]
  
  g_label = ifelse(fiton==0, paste(loci ,"Drift",pop_s, gen_s, sep="_" ), paste(loci, epi, pop_s, gen_s, sep="_"))
  
  for (i in 1:(length(f_list))){
    ## collate al_freq files
    filename = f_list [i]
    al_freq  = fread(file = paste0(path,"group_",g, "/",filename), fill=TRUE,skip="Gen")
    al_freq [,run := i]
    al_freq [,group:=g]
    al_freq [,y:=epi]
    al_freq [,l:=loci]
    al_freq [,fit:=fiton]
    al_freq [,s_gen:=sum_gen]
    al_freq [,w_gen:=win_gen]
    al_freq [,s_pop:=sum_pop]
    al_freq [,w_pop:=win_pop]
    al_freq[, pop_season := pop_s]
    al_freq[, gen_season := gen_s]
    al_freq[,label := g_label]
    alfreq_data  = rbind(alfreq_data , al_freq )
  }
}
alfreq_data[gen_season == "EG", gen_year:=Gen%%20, by=c("label", "Gen")] # set generation in the seasonal cycle/year
alfreq_data[gen_season == "UG", gen_year:=Gen%%10, by=c("label", "Gen")] 
alfreq_data[gen_season == "EG", season:=ifelse((gen_year<11 & gen_year>0), "summer", "winter"), by=c("Gen")] # set season
alfreq_data[gen_season == "UG", season:=ifelse((gen_year<11 & gen_year>0), "summer", "winter"), by=c("Gen")]
alfreq_data [, id :=paste0(group, "_", run, "_",mut_pos)] # give each mutation a unique id
alfreq_data[mut_freq!=1, Freq.bin:="Segregating"] # label alleles that are segregating
alfreq_data[mut_freq==1, Freq.bin:="Fixed_Summer"] # label alleles that are fixed
alfreq_data[Freq.bin=="Segregating",n_seg:=.N, by = c("Gen", "label","l")] # count the number of segregating alleles
alfreq_data[, env:=paste(pop_season, gen_season, sep="_"), by="label"] # set env variable (for if using unequal generations per season and fluctuating population size)
alfreq_data[, year:=ifelse(gen_season=="EG",Gen%/%20, Gen%/%12), by=c("Gen", "gen_season")] # set year/seasonal cycle
alfreq_data[, sampTime:=ifelse(gen_season=="EG",(20*year)%/%9000, (12*year)%/%9000) , by=c("label","year")] # sampling time
alfreq_data[Freq.bin=="Segregating",fc_year:=max(mut_freq)-min(mut_freq), by=c("id", "year")] # yearly frequency change of segregating alleles
alfreq_data[Gen>=10000, fluctuating:=ifelse(mut_freq!=0 & mut_freq != 1, T, F), by="id"] # determine if fluctuating after initial establishment period 

## Figure S1 - Distribution of mean effect size for segregating and fixed/lost alleles.
lostloc=unique(alfreq_data[Gen==18000&Freq.bin=="Segregating", .(id, y, l)]) #videntify loci segregating at final sampling gen
vals=na.omit(alfreq_data[Gen==1, .(Gen, id, s_d, w_d, s_fx, w_fx, y,l)]) # remove NA values
vals[, finalseg:="Fixed"] # give all "Fixed" label
vals[id %in% lostloc$id, finalseg:="Segregating"] # labels segregating loci as segregating
vals[, meanfx:= (s_fx/w_fx), by=c("id","y", "l")] # average effect size
vals[, y_lab:=paste0("y = ",y), by="y"] # add y label
vals[, l_lab:=paste0("l = ",l), by="l"] # add l label

y_labels=c("y = 0.5", "y = 1","y = 2","y = 4","y = 8", "y = 12", "y = 16","y = 20")
names(y_labels)<-c("0.5", "1", "2", "4","8","12","16","20")

plot=ggplot()+ geom_density(data=vals[finalseg=="Fixed"], aes(x=meanfx), col="red")+ 
  geom_density(data=vals[finalseg=="Segregating"], aes(x=meanfx))+ facet_wrap(~l_lab+y_lab) +
  scale_x_log10()+ labs(x="Seasonal Effect Size Ratio", y="Density") + theme_bw()
ggsave(filename =paste0("plots/FigureS1.jpg"), plot = plot , width = 8, height = 6)
ggsave(filename =paste0("plots/FigureS1.pdf"), plot = plot , width = 8, height = 6)

## Figure 1 - The average number of segregating loci across time.
 
rm.col=c("s_d", "w_d", "s_fx", "w_fx") # remove dominance and effect size columns
alfreq_data[, (rm.col):=NULL]
 # calculate minimum and maximum frequency in year/seasonal cycle
intermediate = alfreq_data[Freq.bin=="Segregating", .(ymax = max(mut_freq), ymin = min(mut_freq)), by=c("id", "year","sampTime","group", "label", "env", "y", "l") ]
 # calculate amplitude across year
intermediate[, yamp:=ymax-ymin, by=c("id", "year","sampTime","group", "label") ]
 # calculate max, min and mean amplitude across replicates
amplitudes = intermediate[yamp>0, .(av_max = mean(ymax), av_min = mean(ymin), av_amp = mean(yamp)), by=c("id","group", "label", "sampTime", "env", "y", "l")]
 # average amplitude for each parameter set
amplitudes[, amp := av_max-av_min, c("id", "group","label", "sampTime", "env", "y")]
amplitudes[, year:=as.integer(sampTime)*9000, by="sampTime"]  # set year
amplitudes[, clean_lab:=paste0("l = ", l), by= "l"] # add l label

segloc_data=unique(na.omit(alfreq_data[Gen>1, .(Gen, n_seg, y, id, env, l)])) ## isolate unique mutation information
segloc_data[, n_seg_rep:=n_seg/20, by=c("Gen", "y", "l" )] # avergae the number of segregating loci across replicates
segdata=unique(segloc_data[Gen<28000, .(Gen, y, l, n_seg_rep)]) #remove any additional duplicates
segdata[, clean_l:=paste0("l = ", l), by = "l"] # add l label
# create empty data.table for figure
dummy <- data.table(Gen = c(0, 27000, 0, 27000, 0, 27000 ), n_seg = c(0, 100, 0, 200, 0, 500), l=c(100,100,200,200,500,500))

l_labels<-c("l = 100", "l = 200", "l = 500") # set labels
names(l_labels)<-c("100", "200", "500")
Fig1 = ggplot(segdata, aes(x = Gen)) + 
  geom_blank(data=dummy, aes(x=Gen, y=n_seg))+
  geom_line(aes(y=n_seg_rep,  col=factor(y)), alpha=0.8) +
  theme_bw()+ coord_cartesian(xlim=c(0,17000))+
  labs(x="Generations",y="Number of Segregating Loci", col="Epistasis")+facet_wrap(~l, scale="free_y", labeller = labeller(l=l_labels))

ggsave(filename =paste0("plots/Figure1.pdf"), plot = Fig1 , width = 11, height = 6)
ggsave(filename =paste0("plots/Figure1.jpg"), plot = Fig1 , width = 11, height = 6)

## Figure 2 - Distribution of the amplitude of allele frequency fluctuations across time. 

amplitudes[, clean_y:=paste0("y = ", y), by= "y"] # add y label

al.dist = ggplot(amplitudes[y>1 & year<20000 ],aes(x = amp, col = factor(year))) + 
  geom_density() + scale_x_log10()+
  theme_bw()+
  labs(x="Amplitude of Fluctuations", y="Density", col= "Year")+ facet_wrap(~clean_lab+factor(clean_y, levels = c("y = 2","y = 4", "y = 8", "y = 12", "y = 16", "y = 20")), nrow=6)

ggsave(filename="plots/Figure2.pdf", al.dist, width = 8, height =8)
ggsave(filename="plots/Figure2.jpg", al.dist, width = 8, height =8)

## Figure 3 - Relationship between number of segregating loci, epistasis, and amplitude of seasonal allele frequency fluctuation.
# calculate average max, mean and minimum allele frequency
# calculate minimum and maximum frequency in year/seasonal cycle
intermediate = alfreq_data[Freq.bin=="Segregating", .(ymax = max(mut_freq), ymin = min(mut_freq)), by=c("id", "year","sampTime","group", "label", "env", "y", "l", "run") ]
# calculate amplitude across year
intermediate[, yamp:=ymax-ymin, by=c("id", "year","sampTime","group", "label", "run") ]

amps = intermediate[yamp>0, .(av_max = mean(ymax), av_min = mean(ymin), av_amp = mean(yamp)), by=c("id","group", "label", "sampTime", "env", "y", "l","run")]
# calculate the average amplitude per replicate
amps[, amp := av_max-av_min, c("id", "group","label", "sampTime", "env", "y","run")]
amps[, year:=as.integer(sampTime)*9000, by="sampTime"] # set year
nseg.scaled=amps[, .N, by=c("label","run", "year", "env","l","y")] # count the number of segregating loci at each time point
# calculate the mean, variance, median and 90% quantile of amplitude for each replicate
amps[, `:=` (mean=mean(amp), var=var(amp), med=median(amp), nq=quantile(amp, probs=0.9)), by=c("label","run", "year", "env","y")]
# merge data sets for plotting
nseg.scaled=merge(nseg.scaled, unique(amps[, .(group, label, run, year, mean, var, med, nq)]), by=c("label", "run", "year"))

# rename mean, median and 90% Quartile labels
xx <- names(nseg.scaled); xx[9] <- "Mean"; xx[11] <-"Median"; xx[12] <-"90% Quartile"
setnames(nseg.scaled, xx)
# transform data set to plot in one figure
fig3.data=melt(nseg.scaled[,c(1:9,11:12)], measure.vars = c("Mean", "Median", "90% Quartile"), variable.name = "stat")

fig3.data[, clean_l:=paste0("l = ", l), by= "l"] # set labels
fig3.data[, clean_y:=paste0("y = ", y), by= "y"]

fig.3=ggplot(fig3.data[year==18000],aes(x=N, y=value,col=factor(y)))+
  geom_point(aes(shape=factor(l)))+
  geom_smooth(aes(col=factor(y)),method='lm', formula= y~x, se = FALSE, alpha=0.3)+  
  facet_wrap(~stat, ncol=1, strip.position = 'right')+
  theme_bw()+
  labs(x="Final Number of Segregating Loci", y="Amplitude of Allele Frequency Fluctuation", col="Epistasis (y)", shape=
         "Initial Loci\nNumber")+ scale_x_log10()+ scale_y_log10()

ggsave(filename =paste0("plots/Figure3.jpg"), plot = fig.3 , width = 8, height = 7)
ggsave(filename =paste0("plots/Figure3.pdf"), plot = fig.3 , width = 8, height = 7)

## Multiple regression models and sampling empirically-based simulation parameters

nseg.scaled[, y_numeric := as.numeric(y)] # ensure y values are numeric

# remove the neutral data and select the last year for fitting the model
ID = nseg.scaled$y != "Drift" & nseg.scaled$year == 18000 

# 3D plot, all three axes are logged. Suggests a plain surface would fit well to logged data.
fig <- plot_ly(nseg.scaled[ID ], x = ~N, y = ~Mean, z = ~y_numeric, color = ~y)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'Number Seg Loci', type = 'log'),
                                   yaxis = list(title = 'Mean Amplitude', type = 'log'),
                                   zaxis = list(title = 'Epistasis', type = 'log')))
fig

# Model fit of a multiple regression where all variables are first logged.
nseg.scaled[, log_y := log10(y_numeric)]
nseg.scaled[, log_N := log10(N)]
nseg.scaled[, log_mean := log10(Mean)]

mod1 = lm(log_y ~ log_N + log_mean, data = nseg.scaled[ID])  # Just additive effects
mod2 = lm(log_y ~ log_N*log_mean, data = nseg.scaled[ID])   # With interaction term

summary(mod1)
summary(mod2)

anova(mod1, mod2)  # The model with the interaction term fits significantly better

## Fitting initial loci from final loci
seg.data = nseg.scaled[year==18000] # subset data to just final sampling point
seg.data[, log_l:=log10(l)] # log the loci number
# Plot of 3-dimensions
fig <- plot_ly(seg.data, x = ~l, y = ~N, z = ~y_numeric, color = ~y)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'Inital Seg Loci'),
                                   yaxis = list(title = 'Final Seg Loci'),
                                   zaxis = list(title = 'Epistasis')))
fig
# test both additive and interactive models
smod1= lm(l ~ N+y_numeric, data = seg.data)
smod2= lm(l ~ N*y_numeric, data = seg.data)
summary(smod1)
summary(smod2)
# Test which model fits better, r^2 is slightly higher for interactive model (smod2)
anova(smod1, smod2)

## Derive values from empirical parameters

sample_l=c(111, 90, 264, 27, 9) # final loci numbers from empirical studies (see main text for explanations)
sample_amp=c(0.05, 0.08, 0.02, 0.04, 0.35) # amplitudes of fluctuations from empirical data
samps_mean=data.table(log_N=log10(sample_l), log_mean=log10(sample_amp)) # log the number of final loci and amplitude
mean_y=mod2 %>% predict(samps_mean) # predict the epistasis parameter (y) for each combination of parameters
samps_mean$y=10^mean_y ## un-log y value, used y values rounded to nearest 0.5

samps_loci=data.table(N=sample_l, y_numeric=round(samps_mean$y)) # table of final loci number and derived y values
initial_l=smod2%>% predict(samps_loci) # predict the initial loci number
samps_loci$L=initial_l # add initial loci number to table


