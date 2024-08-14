### Chapter 4 ###


##Collate allele frequencies
library(data.table)
library(ggpubr)
library(viridis)
library(tidyverse)
library(ggdist)
library(readxl)
library(grDevices)
library(cowplot)

path="~/FluctuatingSelectionNe/"
setwd(path)
groups=c(1:x) # parameter set IDs

alfreq_data = NULL
dropcols=c("s_d", "w_d", "s_fx", "w_fx")

for (g in groups){
  print(g)
  
  f_list  <- list.files(path =paste0(path,"/group_", g, "/"),pattern ="al_freq_")
  
  
  parameters <- fread(file=paste0(path,"parameters/group_",g, ".txt"), sep = ":")
  setkey(parameters, V1)
  
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
alfreq_data[gen_season == "EG", gen_year:=Gen%%20, by=c("label", "Gen")]
alfreq_data[gen_season == "UG", gen_year:=Gen%%10, by=c("label", "Gen")] 
alfreq_data[gen_season == "EG", season:=ifelse((gen_year<11 & gen_year>0), "summer", "winter"), by=c("Gen")]
alfreq_data[gen_season == "UG", season:=ifelse((gen_year<11 & gen_year>0), "summer", "winter"), by=c("Gen")]
alfreq_data [, id :=paste0(group, "_", run, "_",mut_pos)]
alfreq_data[mut_freq!=1, Freq.bin:="Segregating"]
alfreq_data[mut_freq==1, Freq.bin:="Fixed_Summer"]
alfreq_data[Freq.bin=="Segregating",n_seg:=.N, by = c("Gen", "label","l")]
sd_val = na.omit(alfreq_data[,.(s_d, w_d, s_fx, w_fx), by=c("id")])
rm.col=c("s_d", "w_d", "s_fx", "w_fx")
alfreq_data[, (rm.col):=NULL]
alfreq_data[, env:=paste(pop_season, gen_season, sep="_"), by="label"]

alfreq_data[, year:=ifelse(gen_season=="EG",Gen%/%20, Gen%/%12), by=c("Gen", "gen_season")]
alfreq_data[, sampTime:=ifelse(gen_season=="EG",(20*year)%/%9000, (12*year)%/%9000) , by=c("label","year")]
alfreq_data[Freq.bin=="Segregating",fc_year:=max(mut_freq)-min(mut_freq), by=c("id", "year")]
alfreq_data[Gen>=10000, fluctuating:=ifelse(mut_freq!=0 & mut_freq != 1, T, F), by="id"]
save(alfreq_data, file=paste0(path,"fluctuating_allele_frequencies.RData"))

load(file=paste0(path,"fluctuating_allele_frequencies.RData"))


intermediate = alfreq_data[Freq.bin=="Segregating", .(ymax = max(mut_freq), ymin = min(mut_freq)), by=c("id", "year","sampTime","group", "label", "env", "y", "l") ]
intermediate[, yamp:=ymax-ymin, by=c("id", "year","sampTime","group", "label") ]
amplitudes = intermediate[yamp>0, .(av_max = mean(ymax), av_min = mean(ymin), av_amp = mean(yamp)), by=c("id","group", "label", "sampTime", "env", "y", "l")]
amplitudes[, amp := av_max-av_min, c("id", "group","label", "sampTime", "env", "y")]
amplitudes[, year:=as.integer(sampTime)*9000, by="sampTime"]
numloci=amplitudes[, .N, by=c("label", "year","env", "y")]
amplitudes[, clean_lab:=paste0("l = ", l), by= "l"]

al.dist = ggplot(amplitudes[ y>1 ],aes(x = as.character(year),y = amp, col = factor(y))) + 
  geom_density(alpha=0,width = 0.8) +
   coord_cartesian(y=c(0,1.0)) +
  theme_bw()+
  scale_x_discrete(limits=c("0","9000","18000"))+
  labs(x="Sampling time", y="Amplitude of Fluctuations", col= "Epistasis")+ facet_wrap(~clean_lab+factor(y, levels = c(2,4, 8, 12, 16, 20)), ncol=6)
ggexport(al.dist, filename=paste0("plots/allele_dist_time_", e.con,"_upscaled.pdf"), width = 15, height =8)
ggsave(filename="plots/fig_1_18kgen.pdf", al.dist, width = 10, height =8)

nloc=unique(alfreq_data$l)
segloc_data=unique(na.omit(alfreq_data[Gen>1, .(Gen, n_seg, y, id, env, l)]))
segloc_data[, n_seg_rep:=n_seg/20, by=c("Gen", "y", "l" )]
segdata=unique(segloc_data[Gen<28000, .(Gen, y, l, n_seg_rep)])
segdata[, clean_l:=paste0("l = ", l), by = "l"]
dummy <- data.table(Gen = c(0, 27000, 0, 27000, 0, 27000 ), n_seg = c(0, 100, 0, 200, 0, 500), l=c(100,100,200,200,500,500))

l_labels<-c("l = 100", "l = 200", "l = 500")
names(l_labels)<-c("100", "200", "500")
Fig1 = ggplot(segdata, aes(x = Gen)) + 
  geom_blank(data=dummy, aes(x=Gen, y=n_seg))+
  geom_line(aes(y=n_seg_rep,  col=factor(y)), alpha=0.8) +
  #ylim(y=c(0,20*nloc)) +
  theme_bw()+ coord_cartesian(xlim=c(0,17000))+
  labs(x="Generations",y="Number of Segregating Loci", col="Epistasis")+facet_wrap(~l, scale="free_y", labeller = labeller(l=l_labels))


# seg.loci
ggexport(seg.loci, filename=paste0("plots/n_seg_L",nloc,"_upscaled.pdf"), width = 15, height =8)
ggsave(filename =paste0("plots/n_seg_L",nloc,"_upscaled.jpg"), plot = seg.loci , width = 11, height = 6)
ggsave(filename =paste0("plots/n_seg_scaled.pdf"), plot = Fig1 , width = 11, height = 6)

amps[, clean_l:=paste0("l = ", l), by= "l"]
amps[, clean_y:=paste0("y = ", y), by= "y"]

e.con="CP_EG"
al.dist = ggplot(amps[env%in%e.con & y>1 & year<20000 ],aes(x = amp, col = factor(year))) + 
  geom_density() + scale_x_log10()+
  theme_bw()+
  labs(x="Amplitude of Fluctuations", y="Density", col= "Year")+ facet_wrap(~clean_lab+factor(clean_y, levels = c("y = 2","y = 4", "y = 8", "y = 12", "y = 16", "y = 20")), nrow=6)
# al.dist
ggexport(al.dist, filename=paste0("plots/allele_dist_time_", e.con,"_upscaled.pdf"), width = 15, height =8)
ggsave(filename="plots/fig2.pdf", al.dist, width = 8, height =8)


lostloc=unique(alfreq_data[Gen==18000&Freq.bin=="Segregating", .(id, y, l)])
# vals=na.omit(alfreq_data[Gen==1 &id %in% lostloc$id, .(s_fx, w_fx, s_d, w_d), by="id"])
freq=alfreq_data[Gen==29, .(fc_year, y, l), by="id"]
vals=na.omit(alfreq_data[Gen==1, .(Gen, id, s_d, w_d, s_fx, w_fx, y,l)])
vals[, finalseg:="Fixed"]
vals[id %in% lostloc$id, finalseg:="Segregating"]
data=merge(vals, freq, by=c("id", "y", "l"))
N.loci=lostloc[, .N, by=c("y","l")]
N.loci[, av:=(N/20)]
data=alfreq_data[Gen<30 & Gen>=20 &id %in% lostloc$id]
data[, fc_year:=max(mut_freq)-min(mut_freq), by=c("id")]
data=data[Gen==29]
rm.col=c("s_d", "w_d", "s_fx", "w_fx")
data[, (rm.col):=NULL]
data=merge(data, vals,by = "id" )
ggplot(data, aes(x=s_fx/w_fx, y=fc_year, col=finalseg))+ geom_point(alpha=0.5)+ facet_wrap(~l+y)+ scale_x_log10()
ggplot(data, aes(x=s_d/w_d, y=fc_year, col=finalseg))+ geom_point(alpha=0.5)+ facet_wrap(~l+y) + scale_x_log10() + labs(x="Mean Dominance", y="Frequency change at the end of first sampling point")


plot=ggplot(data[finalseg=="Fixed"], aes(y=s_d/w_d, x=s_fx/w_fx, col=finalseg))+geom_abline(slope=12)+ geom_point(alpha=0.3)+ facet_wrap(~l.y+y.y) + scale_x_log10()+ scale_y_log10()

data[, meanfx:= (s_fx/w_fx), by=c("id","y", "l")]
data[, y_lab:=paste0("y = ",y), by="y"]
data[, l_lab:=paste0("l = ",l), by="l"]

y_labels=c("y = 0.5", "y = 1","y = 2","y = 4","y = 8", "y = 12", "y = 16","y = 20")
names(y_labels)<-c("0.5", "1", "2", "4","8","12","16","20")

plot=ggplot()+ geom_density(data=data[finalseg=="Fixed"], aes(x=meanfx), col="red")+ 
  geom_density(data=data[finalseg=="Segregating"], aes(x=meanfx))+ facet_wrap(~l_lab+y_lab) +
  scale_x_log10()+ labs(x="Mean Effect Size", y="Density") + theme_bw()
ggsave(filename =paste0("plots/non_seg_fx.jpg"), plot = plot , width = 8, height = 6)



alfreq_data[Gen>=100000, fluctuating:=ifelse(mut_freq!=0 & mut_freq != 1, T, F), by="id"]
seg.alleles=alfreq_data[fluctuating==T&run==5]
e.con="FP_EG"
nloc=unique(alfreq_data$l)

plotloci=c()

for(i in c("Drift", 0.5, 1,2,4, 8, 12, 16, 20)){
  if (length(unique(seg.alleles[env%in%e.con& y==i,id]))<=20){ 
    plotloci=c(plotloci, unique(seg.alleles[env%in%e.con& y==i,id]))
  }else{
  samp.alleles=sample(unique(seg.alleles[env%in%e.con& y==i, id]), 20)
  plotloci=c(plotloci, samp.alleles)}
}

al.freq =ggplot(seg.alleles[id%in%plotloci &env%in%e.con&Gen>=99000 &Gen<100000], aes(x = Gen)) + 
  geom_line(aes(y=mut_freq, group=id, col=id), alpha=0.3) +
  coord_cartesian(y=c(0,1.0)) +
  theme_bw()+
  # scale_x_discrete(limits=c("0","9000","18000","27000", "36000", "45000", "54000", "63000", "72000", "81000", "90000", "99000"), breaks=c("0","9000","18000","27000", "36000", "45000", "54000", "63000", "72000", "81000", "90000", "99000"))+
  ylab("Allele Frequency")+
  # geom_text(data = numloci, aes(x = as.character(year), label=paste0("n=",N), y=.95), size = 4)+
  xlab("Sampling time")+ facet_wrap(env~factor(y, levels = c("Drift",0.5, 1, 2, 4, 8, 12, 16, 20)), ncol=4)+theme(legend.position = 'None')
ggexport(al.freq, filename=paste0("plots/af1r_", nloc,"_",e.con, ".pdf"), width = 15, height =8)


numloci=amplitudes[(data_type=="Simulation" & year==18000)|(data_type=="Empirical"), .N, by=c("label", "year","env", "y", "l")]
numloci[, clean_l:=paste0("l = ", l), by= "l"]
numloci[, clean_y:=paste0("y = ", y), by= "y"]

ampfit = ggplot(amplitudes, aes(x = y)) + 
  geom_jitter(aes(y = amp,col = factor(y)), alpha=0.2, width = 0.38)  +
  geom_boxplot(aes(y = amp),alpha=0,width = 0.8, outlier.color = NA) +
  coord_cartesian(y=c(0,1.0)) +
  theme_bw()+
  ylab("Amplitude of Fluctuations")+
  geom_text(data = numloci, aes(x =y,label=paste0("n=",N), y=-0.025), size = 2)+
   scale_x_discrete(limits=c( 0.5, 1,2,4,8,12,16,20,"Bergland", "Machado"),
                    breaks=c( "0.5", "1","2","4","8","12","16","20","Bergland", "Machado"))+
  xlab("Epistasis")+ facet_wrap(clean_l~env)+theme(legend.position = "None")

ggexport(ampfit, filename=paste0("plots/amp_fit_l",nloc,".pdf"), width = 12, height =6)
ggsave(filename =paste0("plots/amp_fit_", nloc, ".jpg"), plot = ampfit , width = 11, height = 6)



## slope
setwd("~/phd_data/Results/chp4/")

load(file=paste0("~/Box Sync/data/ch4_allele_l1000.RData"))


intermediate = alfreq_data[Freq.bin=="Segregating", .(ymax = max(mut_freq), ymin = min(mut_freq)), by=c("id", "year","sampTime","group", "label", "env", "y", "run") ]
intermediate[, yamp:=ymax-ymin, by=c("id", "year","sampTime","group", "label", "run") ]
amps = intermediate[yamp>0, .(av_max = mean(ymax), av_min = mean(ymin), av_amp = mean(yamp)), by=c("id","group", "label", "sampTime", "env", "y", "run")]
amps[, amp := av_max-av_min, c("id", "group","label", "sampTime", "env", "y", "run")]
amps[, year:=as.integer(sampTime)*9000, by="sampTime"]
nseg.1000=amps[, .N, by=c("label","run", "year", "env","y")]
amps[, `:=` (mean=mean(amp), var=var(amp), med=median(amp), nq=quantile(amp, probs=0.9)), by=c("label","run", "year", "env","y")]
nseg.1000=merge(nseg.1000, unique(amps[, .(group, label, run, year, mean, var, med, nq)]), by=c("label", "run", "year"))
nseg.1000[, loci:=1000]

load(file=paste0("~/Box Sync/data/ch4_allele_l500.RData"))


intermediate = alfreq_data[Freq.bin=="Segregating", .(ymax = max(mut_freq), ymin = min(mut_freq)), by=c("id", "year","sampTime","group", "label", "env", "y", "run") ]
intermediate[, yamp:=ymax-ymin, by=c("id", "year","sampTime","group", "label", "run") ]
amps = intermediate[yamp>0, .(av_max = mean(ymax), av_min = mean(ymin), av_amp = mean(yamp)), by=c("id","group", "label", "sampTime", "env", "y", "run")]
amps[, amp := av_max-av_min, c("id", "group","label", "sampTime", "env", "y", "run")]
amps[, year:=as.integer(sampTime)*9000, by="sampTime"]
nseg.500=amps[, .N, by=c("label","run", "year", "env","y")]
amps[, `:=` (mean=mean(amp), var=var(amp), med=median(amp), nq=quantile(amp, probs=0.9)), by=c("label","run", "year", "env","y")]
nseg.500=merge(nseg.500, unique(amps[, .(group, label, run, year, mean, var, med, nq)]), by=c("label", "run", "year"))
nseg.500[, loci:=500]

load(file=paste0("~/Box Sync/data/ch4_allele_l100.RData"))


intermediate = alfreq_data[Freq.bin=="Segregating", .(ymax = max(mut_freq), ymin = min(mut_freq)), by=c("id", "year","sampTime","group", "label", "env", "y", "run") ]
intermediate[, yamp:=ymax-ymin, by=c("id", "year","sampTime","group", "label", "run") ]
amps = intermediate[yamp>0, .(av_max = mean(ymax), av_min = mean(ymin), av_amp = mean(yamp)), by=c("id","group", "label", "sampTime", "env", "y", "run")]
amps[, amp := av_max-av_min, c("id", "group","label", "sampTime", "env", "y", "run")]
amps[, year:=as.integer(sampTime)*9000, by="sampTime"]
nseg.100=amps[, .N, by=c("label","run", "year", "env","y")]
amps[, `:=` (mean=mean(amp), var=var(amp), med=median(amp), nq=quantile(amp, probs=0.9)), by=c("label","run", "year", "env","y")]
nseg.100=merge(nseg.100, unique(amps[, .(group, label, run, year, mean, var, med, nq)]), by=c("label", "run", "year"))
nseg.100[, loci:=100]

load(file=paste0("~/Box Sync/data/ch4_allele_l200.RData"))


intermediate = alfreq_data[Freq.bin=="Segregating", .(ymax = max(mut_freq), ymin = min(mut_freq)), by=c("id", "year","sampTime","group", "label", "env", "y", "run") ]
intermediate[, yamp:=ymax-ymin, by=c("id", "year","sampTime","group", "label", "run") ]
amps = intermediate[yamp>0, .(av_max = mean(ymax), av_min = mean(ymin), av_amp = mean(yamp)), by=c("id","group", "label", "sampTime", "env", "y", "run")]
amps[, amp := av_max-av_min, c("id", "group","label", "sampTime", "env", "y", "run")]
amps[, year:=as.integer(sampTime)*9000, by="sampTime"]
nseg.200=amps[, .N, by=c("label","run", "year", "env","y")]
amps[, `:=` (mean=mean(amp), var=var(amp), med=median(amp), nq=quantile(amp, probs=0.9)), by=c("label","run", "year", "env","y")]
nseg.200=merge(nseg.200, unique(amps[, .(group, label, run, year, mean, var, med, nq)]), by=c("label", "run", "year"))
nseg.200[, loci:=200]

# nseg.all=rbind(nseg.100,  nseg.200,nseg.500, nseg.1000)


envs="CP_EG"
nseg.all=rbind(nseg.100[env%in%envs],  nseg.200[env%in%envs],nseg.500[env%in%envs])#, nseg.1000[env%in%envs])

mean.amp=ggplot(nseg.all[year==99000],aes(x=N, y=mean))+
  geom_point(aes(col=factor(loci)))+
  geom_smooth(method='loess', formula= y~x, se = FALSE)+
  theme_bw()+
  labs(x="Final N segregating", y="Mean amplitude")+ 
  # geom_violin()+
  facet_wrap(env~factor(y, levels = c("Drift",0.5, 1, 2, 4, 8, 12, 16, 20)), scales = "free")
ggexport(mean.amp, filename=paste0("plots/Scaling_meanamp_",envs,".pdf"), width = 12, height =6)

var.amp=ggplot(nseg.all[year==99000],aes(x=N, y=var))+
  geom_point(aes(col=factor(loci)))+
  geom_smooth(method='loess', formula= y~x, se = FALSE)+
  theme_bw()+
  labs(x="Final N segregating", y="Variance in amplitude")+
  # geom_violin()+
  facet_wrap(env~factor(y, levels = c("Drift",0.5, 1, 2, 4, 8, 12, 16, 20)), scales = "free")
ggexport(var.amp, filename=paste0("plots/Scaling_varamp_",envs,".pdf"), width = 12, height =6)

med.amp=ggplot(nseg.all[year==99000],aes(x=N, y=med))+
  geom_point(aes(col=factor(loci)))+
  geom_smooth(method='loess', formula= y~x, se = FALSE)+
  theme_bw()+
  labs(x="Final N segregating", y="Median amplitude")+ 
  # geom_violin()+
  facet_wrap(env~factor(y, levels = c("Drift",0.5, 1, 2, 4, 8, 12, 16, 20)), scales = "free")
ggexport(med.amp, filename=paste0("plots/Scaling_medamp_",envs,".pdf"), width = 12, height =6)


slope=ggplot(nseg.all[year==99000],aes(x=N, y=med))+
  geom_point(aes(shape=factor(loci),col=factor(y)))+
  geom_smooth(aes(col=factor(y)),method='lm', formula= y~x, se = FALSE)+
  theme_bw()+
  labs(x="Final N segregating", y="Median amplitude")+ scale_x_log10()+ scale_y_log10()
ggexport(slope, filename=paste0("plots/slope_med_",envs,".pdf"), width = 12, height =6)

slope=ggplot(nseg.all[year==99000],aes(x=N, y=mean))+
  geom_point(aes(shape=factor(loci),col=factor(y)))+
  geom_smooth(aes(col=factor(y)),method='lm', formula= y~x, se = FALSE)+
  theme_bw()+
  labs(x="Final N segregating", y="Mean amplitude")+ scale_x_log10()+ scale_y_log10()
ggexport(slope, filename=paste0("plots/slope_mean_",envs,".pdf"), width = 12, height =6)

slope=ggplot(nseg.all[year==99000],aes(x=N, y=nq))+
  geom_point(aes(shape=factor(loci),col=factor(y)))+
  geom_smooth(aes(col=factor(y)),method='lm', formula= y~x, se = FALSE)+
  theme_bw()+
  labs(x="Final N segregating", y="90th percentile amplitude")+ scale_x_log10()+ scale_y_log10()
ggexport(slope, filename=paste0("plots/slope_nq_",envs,".pdf"), width = 12, height =6)


## SCALED COMPARISONS


load(file=paste0("~/Box Sync/data/ch4_allele_scaled.RData"))

intermediate = alfreq_data[Freq.bin=="Segregating", .(ymax = max(mut_freq), ymin = min(mut_freq)), by=c("id", "year","sampTime","group", "label", "l", "env", "y", "run") ]
intermediate[, yamp:=ymax-ymin, by=c("id", "year","sampTime","group", "label","run") ]
amps = intermediate[yamp>0, .(av_max = mean(ymax), av_min = mean(ymin), av_amp = mean(yamp)), by=c("id","group", "label", "sampTime", "env", "y", "l","run")]
amps[, amp := av_max-av_min, c("id", "group","label", "sampTime", "env", "y","run")]
amps[, year:=as.integer(sampTime)*9000, by="sampTime"]
nseg.scaled=amps[, .N, by=c("label","run", "year", "env","l","y")]
amps[, `:=` (mean=mean(amp), var=var(amp), med=median(amp), nq=quantile(amp, probs=0.9)), by=c("label","run", "year", "env","y")]
nseg.scaled=merge(nseg.scaled, unique(amps[, .(group, label, run, year, mean, var, med, nq)]), by=c("label", "run", "year"))

for ( i in c(1,20)){
  scaled=nseg.scaled[y==i & year==27000]
  init=nseg.100[y==i& year==54000 & env=="CP_EG"]
  t.test(scaled[,mean], init[,mean])
  t.test(scaled[,med], init[,med])
  t.test(scaled[,nq], init[,nq])
}

median.slope=ggplot(nseg.scaled[year==18000],aes(x=N, y=med,col=factor(y)))+
  geom_point(aes(shape=factor(l)))+
  geom_smooth(aes(col=factor(y)),method='lm', formula= y~x, se = FALSE, alpha=0.5)+
  theme_bw()+
  labs(x="Final N segregating", y="Median amplitude")+ scale_x_log10()+ scale_y_log10()+theme(axis.title.x = element_blank())
slope
ggexport(slope, filename=paste0("plots/slope_med_scaled_log.pdf"), width = 12, height =6)

mean.slope=ggplot(nseg.scaled[year==18000],aes(x=N, y=mean,col=factor(y)))+
  geom_point(aes(shape=factor(l)))+
  geom_smooth(aes(col=factor(y)),method='lm', formula= y~x, se = FALSE, alpha=0.5)+
  theme_bw()+
  labs(x="Final N segregating", y="Mean amplitude")+ scale_x_log10()+ scale_y_log10()+theme(axis.title.x = element_blank())
ggexport(slope, filename=paste0("plots/slope_mean_scaled_nolog.pdf"), width = 12, height =6)

nq.slope=ggplot(nseg.scaled[year==18000],aes(x=N, y=nq,col=factor(y)))+
  geom_point(aes(shape=factor(l)))+
  geom_smooth(aes(col=factor(y)),method='lm', formula= y~x, se = FALSE, alpha=0.3)+  
  # geom_point(data=Q9, aes(y=amp,col=y))+
  theme_bw()+
  labs(x="Final Number of Segregating Loci", y="90th percentile")+ scale_x_log10()+ scale_y_log10()+ theme(legend.position='none')
ggexport(slope, filename=paste0("plots/slope_nq_scaled_nolog.pdf"), width = 12, height =6)

xx <- names(nseg.scaled); xx[9] <- "Mean"; xx[11] <-"Median"; xx[12] <-"90% Quartile"
setnames(nseg.scaled, xx)

fig3.data=melt(nseg.scaled[,c(1:9,11:12)], measure.vars = c("Mean", "Median", "90% Quartile"), variable.name = "stat")

fig3.data[, clean_l:=paste0("l = ", l), by= "l"]
fig3.data[, clean_y:=paste0("y = ", y), by= "y"]

fig.3=ggplot(fig3.data[year==18000],aes(x=N, y=value,col=factor(y)))+
  geom_point(aes(shape=factor(l)))+
  geom_smooth(aes(col=factor(y)),method='lm', formula= y~x, se = FALSE, alpha=0.3)+  
  facet_wrap(~stat, ncol=1, strip.position = 'right')+
  theme_bw()+
  labs(x="Final Number of Segregating Loci", y="Amplitude of Allele Frequency Fluctuation", col="Epistasis (y)", shape=
         "Initial Loci\nNumber")+ scale_x_log10()+ scale_y_log10()

ggsave(filename =paste0("plots/fig3.jpg"), plot = fig.3 , width = 8, height = 7)


fig3.data[year==18000, .(min=range(value)[1], max=range(value)[2]), by=c("stat", "l", "y")]

