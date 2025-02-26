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
amps[, `:=` (mean=mean(amp), var=var(amp), med=median(amp), max=max(amp),nq=quantile(amp, probs=0.9)), by=c("label","run", "year", "env","y")]
# merge data sets for plotting
nseg.scaled=merge(nseg.scaled, unique(amps[, .(group, label, run, year, mean, var, med,max, nq)]), by=c("label", "run", "year"))

# rename mean, median,max  and 90% Quartile labels
xx <- names(nseg.scaled); xx[9] <- "Mean"; xx[11] <-"Median"; xx[12] <-"90% Quartile"; xx[13] <- "Maximum"
setnames(nseg.scaled, xx)
# transform data set to plot in one figure
fig1.data=melt(nseg.scaled[,c(1:9,11:13)], measure.vars = c("Mean", "Median", "Maximum","90% Quartile"), variable.name = "stat")

fig1.data[, clean_l:=paste0("l = ", l), by= "l"] # set labels
fig1.data[, clean_y:=paste0("y = ", y), by= "y"]

fig.1=ggplot(fig1.data[year==18000],aes(x=N, y=value,col=factor(y)))+
  geom_point(aes(shape=factor(l)))+
  geom_smooth(aes(col=factor(y)),method='lm', formula= y~x, se = FALSE, alpha=0.3)+  
  facet_wrap(~stat, ncol=1, strip.position = 'right')+
  theme_bw()+
  labs(x="Final Number of Segregating Loci", y="Amplitude of Allele Frequency Fluctuation", col="Epistasis (y)", shape=
         "Initial Loci\nNumber")+ scale_x_log10()+ scale_y_log10()

ggsave(filename =paste0("plots/Figure1.pdf"), plot = fig.1 , width = 11, height = 6)
ggsave(filename =paste0("plots/Figure1.jpg"), plot = fig.1 , width = 11, height = 6)

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


