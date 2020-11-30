library(ape)




##Read in probeCertain ASTRAL tree
setwd("~/PleuroceridaePhylogenetics/probeRegion/CertainPruned/multipleIQtreeRuns_finalTrees/")
tree2<-read.tree(file = "probeCertain_LB_RCFV_newIQTREE-BS10.astral.renamed.tre")
tree2$edge.length

##Branch lengths not set by ASTRAL are read in as NaN, which allows use of below command
tree2$edge.length[is.nan(tree2$edge.length)]<-0.15
tree2$edge.length
write.tree(tree2,file="probeCertain_LB_RCFV_newIQTREE-BS10.astral.renamed.rescaled.tre")


##Read in probeCertain astralMap Trees
astralMapTree<-read.tree(file ="probeCertain_LB_RCFV_newIQTREE-BS10.astral.mapCahabaPleuro2Lepto.tre")
astralMapTree$edge.length
astralMapTree$edge.length[is.nan(astralMapTree$edge.length)]<-0.15
write.tree(astralMapTree,file="probeCertain_LB_RCFV_newIQTREE-BS10.astral.mapCahabaPleuro2Lepto.rescaled.tre")

##NEXT TREE  ##Just without taxon map
setwd("~/PleuroceridaePhylogenetics/threeRegions/Alicut/CertainTrimmed/ASTRAL/")
##Read in probeCertain ASTRAL tree
tree2<-read.tree(file = "threeAlicutCertain_LB_RCFV_newIQTREE-BS10.astral.renamed.tre")
tree2$edge.length

##Branch lengths not set by ASTRAL are read in as NaN, which allows use of below command
tree2$edge.length[is.nan(tree2$edge.length)]<-0.15
tree2$edge.length
write.tree(tree2,file="threeAlicutCertain_LB_RCFV_newIQTREE-BS10.astral.renamed.rescaled.tre")

##NEXT TREE
setwd("~/PleuroceridaePhylogenetics/threeRegions/noMasking/CertainTrimmed/multipleIQTREEruns_finalTrees/")
##Read in probeCertain ASTRAL tree
tree2<-read.tree(file = "threeCertain_LB_RCFV_newIQTREE-BS10.astral.renamed.tre")
tree2$edge.length

##Branch lengths not set by ASTRAL are read in as NaN, which allows use of below command
tree2$edge.length[is.nan(tree2$edge.length)]<-0.15
tree2$edge.length
write.tree(tree2,file="threeCertain_LB_RCFV_newIQTREE-BS10.astral.renamed.rescaled.tre")
