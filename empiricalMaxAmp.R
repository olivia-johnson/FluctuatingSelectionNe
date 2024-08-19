## Load packages and install any missing from the user's R library
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table)

## Bergland and Machado Maximum amplitude values

path="~/FluctuatingSelectionNe/"
setwd(path) ## set path to empirical data

## Command to extract desired columns from Bergland et al. dataset (Chrom, Pos, Info, PA time series data)
cmd="cut -f 1,2,8,17,18,19,20,21,22,23 6d_v7.3_output.vcf | sed 's/[\t]/,/g' > bergland.csv"
system(cmd) ## if this gives an error, run the command directly in the command line

## load in Bergland et al. 2014 data (https://datadryad.org/stash/dataset/doi:10.5061/dryad.v883p)
df=fread("bergland.csv") ## if this gives an error, run cmd directly in the command line
## extract precalculated mean spring frequency 
df[,springF:=grep('SprF=', strsplit(INFO, ';')[[1]], value = TRUE), by=c("#CHROM", "POS")]
## extract precalculated mean fall frequency 
df[,fallF:=grep('FallF=', strsplit(INFO, ';')[[1]], value = TRUE), by=c("#CHROM", "POS")]
## extract false discovery rate values
df[,SQ:=grep('SQ=', strsplit(INFO, ';')[[1]], value = TRUE), by=c("#CHROM", "POS")]
## remove info column to speed up analysis
df[,INFO:=NULL]
## transform values to numeric
df[,SQ:=as.numeric(strsplit(SQ, '=')[[1]][2]), by=c("#CHROM", "POS")]
## filter for FDR < 0.3 (as in Bergland et al)
df=df[SQ<0.3]
## transform values to numeric
df[,springF:=as.numeric(strsplit(springF, '=')[[1]][2]), by=c("#CHROM", "POS")]
df[,fallF:=as.numeric(strsplit(fallF, '=')[[1]][2]), by=c("#CHROM", "POS")]
## calculate mean seasonal amplitude
df[,amp:=abs(as.numeric(springF)-as.numeric(fallF)), by=c("#CHROM", "POS")]
## check spring and fall mean values against precalculated values
df[, `:=` (S_09= as.numeric(strsplit(PA_7_2009[1], ":")[[1]][1])/as.numeric(strsplit(PA_7_2009[1], ":")[[1]][2]), 
           S_10= as.numeric(strsplit(PA_7_2010[1], ":")[[1]][1])/as.numeric(strsplit(PA_7_2010[1], ":")[[1]][2]),
           S_11= as.numeric(strsplit(PA_7_2011[1], ":")[[1]][1])/as.numeric(strsplit(PA_7_2011[1], ":")[[1]][2]), 
           F_09= as.numeric(strsplit(PA_11_2009[1], ":")[[1]][1])/as.numeric(strsplit(PA_11_2009[1], ":")[[1]][2]), 
           F_10= as.numeric(strsplit(PA_11_2010[1], ":")[[1]][1])/as.numeric(strsplit(PA_11_2010[1], ":")[[1]][2]), 
           F_11= as.numeric(strsplit(PA_10_2011[1], ":")[[1]][1])/as.numeric(strsplit(PA_10_2011[1], ":")[[1]][2]) 
), by=c("#CHROM", "POS")]

## filter data by alleles that change direction between seasons
df=df[(S_09<F_09 & S_10<F_09 &S_10<F_10 & S_11<F_10)|(S_09>F_09 & S_10>F_09 &S_10>F_10 & S_11>F_10)] 

df[, `:=` (Spr= (S_09+ S_10+S_11)/3, Fall=(F_09+F_10+F_11)/3), by=c("#CHROM", "POS")]
df[springF==signif(Spr,6)]
df[,calc_amp:=abs(Spr-Fall), by=c("#CHROM", "POS")]
## Maximum amplitude for Bergland et al. (0.385)
df[,max(amp)]
df[,max(calc_amp)]

## Load Machado et al. 2021 Rdata object (https://datadryad.org/stash/dataset/doi:10.5061/dryad.4r7b826)
load("mel_freqdp_042016_Ne_fixed_correctBAVI.Rdata")
machado_data=as.data.table(cbind(info, freq)) ## merge chromosome and position information with frequencies
## read in supplementary metadata (https://cdn.elifesciences.org/articles/67577/elife-67577-supp1-v2.xlsx)
m_snps<- as.data.table(read_excel("machado/elife-67577-supp1-v2.xlsx", 
                                  sheet = "Supplementaryfile1B"))
m_pops <-as.data.table(read_excel("machado/elife-67577-supp1-v2.xlsx", 
                                  sheet = "Supplementaryfile1A"))
## Keep Core20 populations used for seasonal analysis
pop_keep = m_pops[Core20=="yes", .(.N, Sample, InternalName, Year, Season), by=c("Locality")]
pops = unique(pop_keep$Locality) ## unique population labels
## Add value for the year of sampling for each locality
for (i in pops){
  pop_keep[Locality==i & min(Year), yr:=ifelse(Year==min(Year), 1, 2)]
}
## labels fall or spring of each
pop_keep[, label:=paste(Locality, yr, Season, sep="_"), by= c("Locality", "yr", "Season")]
## subset Machado data by samples form filters above
mpops= m_pops[Sample%in%pop_keep$Sample]
## change columns to label created earlier
cols=c("X.CHROM", "POS", pop_keep$InternalName)
col.names=c("chrom", "pos", pop_keep$label)
machado=machado_data[, ..cols]
setnames(machado, cols, col.names)
## SNP coordinates 
snp_keep = m_snps[, .(chrom, pos)]
## Merge datasets
machado=merge(snp_keep, machado, all.x = TRUE, by=c("chrom", "pos"))
## make SNP ID
machado[, id:=paste(chrom,pos, sep="_"), by=c("chrom", "pos")]
machado=na.omit(machado) ## remove NA data
## transform data to have each column as a seasonal time point
machado.af = melt(machado, variable.name = "locality", measure.vars = patterns("1_spring", "1_fall", "2_spring","2_fall"), value.name=c("s_1", "f_1", "s_2","f_2"), variable.factor = FALSE)
machado.af=na.omit(machado.af) ## remove NA data
## filter data by alleles that change direction between seasons
macado.af = machado.af[(s_1<f_1 & s_2<f_1 & s_2<f_2)|(s_1>f_1 & s_2>f_1 & s_2>f_2)]
## calculate mean fall (FF) and spring (SF) frequencies to align with Bergland et al. calculations
machado.af[, FF:=((f_1+ f_2)/2), by=c("id", "locality")]
machado.af[, SF:=((s_1+ s_2)/2), by=c("id", "locality")]
## calculate mean amplitude of fluctuations
machado.af[, amp:=abs(FF-SF), by=c("id", "locality")]
## Maximum amplitude of allele frequency from Machado et al. 2021
machado.af[, max(amp)]

## Mean frequency of top 9 SNPs
setorder(machado.af, -amp) ## Machado - order by amplitude
machado.af[1:9, mean(amp)] ## Machado - mean amplitude of top 9
setorder(df, -amp) ## Bergland - order by amplitude
df[1:9, mean(amp)] ## Bergland - mean amplitude of top 9

