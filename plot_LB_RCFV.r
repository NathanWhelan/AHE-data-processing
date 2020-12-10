library(ggplot2)
setwd("~/PleuroceridaePhylogenetics/threeRegions/Alicut/UncertainTrimmed/")  #Set this as appropriate
##This script is for plotting outputs from TreSpEx and BaCoCa. 
##Do no uncritically use the cutoff values used in the below plots, as they're only relevant to the Whelan et al. pleurocerid AHE study.

##############################################
##############LB SCORES#######################
##############################################

##READ in Data for Three Region, Alicut, UncertainPruned
LB_scores<-read.table("threeAlicutUncertain_LB_scores_summary_perPartition.txt", header=TRUE)

##Generate density plots
P_hetero<-ggplot(data=LB_scores, aes(x=LB_score_Heterogeneity)) +geom_density()
P_quartile<-ggplot(data=LB_scores, aes(x=LB_score_upper_quartile))+geom_density()

##plot Graphs with vertical line at proposed cutoff
pdf(file = "threeAlicutUncertain_LB_score_heterogeneity.pdf")
P_hetero+geom_vline(xintercept = 93)+theme_classic() + labs(title="LB Heterogeneity 93 Cutoff", x = "LB Heterogeneity", y= "Density")
dev.off()

pdf(file = "threeAlicutUncertain_LB_score_quartile.pdf")
P_quartile+geom_vline(xintercept = 78)+theme_classic() + labs(title="LB Quartile 78 Cutoff", x = "LB Quartile", y= "Density")
dev.off()


##READ in Data for Three Region, Alicut, CertainPruned
rm(list=ls())
setwd("~/PleuroceridaePhylogenetics/threeRegions/Alicut/CertainTrimmed/")
LB_scores<-read.table("threeAlicutCertain_LB_scores_summary_perPartition.txt", header=TRUE)

##Generate Density Plots
P_hetero<-ggplot(data=LB_scores, aes(x=LB_score_Heterogeneity)) +geom_density()
P_quartile<-ggplot(data=LB_scores, aes(x=LB_score_upper_quartile))+geom_density()

pdf(file = "threeAlicutCertain_LB_score_heterogeneity.pdf")
P_hetero+geom_vline(xintercept = 92)+theme_classic() + labs(title="LB Heterogeneity 92 Cutoff", x = "LB Heterogeneity", y= "Density")
dev.off()

pdf(file = "threeAlicutCertain_LB_score_quartile.pdf")
P_quartile+geom_vline(xintercept = 79)+theme_classic() + labs(title="LB Quartile 79 Cutoff", x = "LB Quartile", y= "Density")
dev.off()


##READ in Data for Three Region, CertainPruned
rm(list=ls())
setwd("~/PleuroceridaePhylogenetics/threeRegions/noMasking/CertainTrimmed/")
LB_scores<-read.table("LB_scores_summary_perPartition.txt", header=TRUE)

##Generate Density Plots
P_hetero<-ggplot(data=LB_scores, aes(x=LB_score_Heterogeneity)) +geom_density()
P_quartile<-ggplot(data=LB_scores, aes(x=LB_score_upper_quartile))+geom_density()

pdf(file = "threeCertain_LB_score_heterogeneity.pdf")
P_hetero+geom_vline(xintercept = 92)+theme_classic() + labs(title="LB Heterogeneity 92 Cutoff", x = "LB Heterogeneity", y= "Density")
dev.off()

pdf(file = "threeCertain_LB_score_quartile.pdf")
P_quartile+geom_vline(xintercept = 90)+theme_classic() + labs(title="LB Quartile 90 Cutoff", x = "LB Quartile", y= "Density")
dev.off()



##READ in Data for ProbeRegion CertainPruned
rm(list=ls())
setwd("~/PleuroceridaePhylogenetics/probeRegion/CertainPruned/")
LB_scores<-read.table("probeRegionCertain_LB_scores_summary_perPartition.txt", header=TRUE)
##Generate Density Plots
P_hetero<-ggplot(data=LB_scores, aes(x=LB_score_Heterogeneity)) +geom_density()
P_quartile<-ggplot(data=LB_scores, aes(x=LB_score_upper_quartile))+geom_density()
pdf(file = "probeRegionCertain_LB_score_heterogeneity.pdf")
P_hetero+geom_vline(xintercept = 96)+theme_classic() + labs(title="LB Heterogeneity 96 Cutoff", x = "LB Heterogeneity", y= "Density")
dev.off()
pdf(file = "probeRegionCertain_LB_score_quartile.pdf")
P_quartile+geom_vline(xintercept = 93)+theme_classic() + labs(title="LB Quartile 93 Cutoff", x = "LB Quartile", y= "Density")
dev.off()

##READ in Data for ProbeRegion UncertainPruned
rm(list=ls())
setwd("~/PleuroceridaePhylogenetics/probeRegion/UncertainPruned/")
LB_scores<-read.table("probeUncertain_LB_scores_summary_perPartition.txt", header=TRUE)
##Generate Density Plots
P_hetero<-ggplot(data=LB_scores, aes(x=LB_score_Heterogeneity)) +geom_density()
P_quartile<-ggplot(data=LB_scores, aes(x=LB_score_upper_quartile))+geom_density()
pdf(file = "probeRegionUnertain_LB_score_heterogeneity.pdf")
P_hetero+geom_vline(xintercept = 81)+theme_classic() + labs(title="LB Heterogeneity 81 Cutoff", x = "LB Heterogeneity", y= "Density")
dev.off()
pdf(file = "probeRegionUnertain_LB_score_quartile.pdf")
P_quartile+geom_vline(xintercept = 95)+theme_classic() + labs(title="LB Quartile 95 Cutoff", x = "LB Quartile", y= "Density")
dev.off()


####################################################################
###############         RCFV         ###############################
####################################################################

rm(list=ls())
##READ in Data for Three Region, Alicut, UncertainPruned
setwd("~/PleuroceridaePhylogenetics/threeRegions/Alicut/UncertainTrimmed/")
RCFV<-read.table("threeAlicutUncertain_summarized_frequencies.txt", header=TRUE, sep = "\t")
P_RCFV<-ggplot(data=RCFV, aes(x=RCFV$RCFV.Value)) +geom_density()
pdf(file="threeAlicutUncertain_RCFV.pdf")
P_RCFV +geom_vline(xintercept=0.0325)+theme_classic() + labs(title="RCVF cutoff 0.0325", x = "RCFV", y= "Density")
dev.off()


##READ in Data for Three Region, Alicut, CertainPruned
rm(list=ls())
setwd("~/PleuroceridaePhylogenetics/threeRegions/Alicut/CertainTrimmed/")
RCFV<-read.table("summarized_frequencies.txt", header=TRUE, sep = "\t")
P_RCFV<-ggplot(data=RCFV, aes(x=RCFV$RCFV.Value)) +geom_density()
pdf(file="threeAlicutCertain_RCFV.pdf")
P_RCFV +geom_vline(xintercept=0.032)+theme_classic() + labs(title="RCVF cutoff 0.032", x = "RCFV", y= "Density")
dev.off()


##READ in Data for Three Region, CertainPruned
rm(list=ls())
setwd("~/PleuroceridaePhylogenetics/threeRegions/noMasking/CertainTrimmed/")
RCFV<-read.table("summarized_frequencies.txt", header=TRUE, sep = "\t")
P_RCFV<-ggplot(data=RCFV, aes(x=RCFV$RCFV.Value)) +geom_density()
pdf(file="threeCertain_RCFV.pdf")
P_RCFV +geom_vline(xintercept=0.037)+theme_classic() + labs(title="RCVF cutoff 0.037", x = "RCFV", y= "Density")
dev.off()





##READ in Data for Probe Region, CertainPruned
setwd("~/PleuroceridaePhylogenetics/probeRegion/CertainPruned/")
rm(list=ls())
RCFV<-read.table("probeRegionCertain_summarized_frequencies.txt", header=TRUE, sep = "\t")
P_RCFV<-ggplot(data=RCFV, aes(x=RCFV$RCFV.Value)) +geom_density()
pdf(file="probeRegionCertain_RCFV.pdf")
P_RCFV +geom_vline(xintercept=0.0225)+theme_classic() + labs(title="RCVF cutoff 0.0225", x = "RCFV", y= "Density")
dev.off()

##READ in Data for Probe Region, UncertainPruned
rm(list=ls())
setwd("~/PleuroceridaePhylogenetics/probeRegion/UncertainPruned/")
RCFV<-read.table("probeUncertain_summarized_frequencies.txt", header=TRUE, sep = "\t")
P_RCFV<-ggplot(data=RCFV, aes(x=RCFV$RCFV.Value)) +geom_density()
pdf(file="probeRegionUncertain_RCFV.pdf")
P_RCFV +geom_vline(xintercept=0.0225)+theme_classic() + labs(title="RCVF cutoff 0.0225", x = "RCFV", y= "Density")
dev.off()
