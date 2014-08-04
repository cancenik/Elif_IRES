library(edgeR)
library (gplots)
library(apcluster)

dat = read.csv ('~/elif_ires/Elif_DataFiles/072314_Elif_comparison_ALLDATA.csv')
dat[,-1] = log10(dat[,-1])

# test whether mean ratio is the across all cell lines
# We can use either kruskal-wallis non-parametric or aov for parametric assumption
cell_type = c("ESC", "ESC", "EB", "EB","EB","EB", "NSC","NSC", "Neuron","Neuron", "Limb", "Limb" ,"Limb" ,"Limb")
kruskal_pvals = apply ( dat[,-1], 1, function (x) {kruskal.test (as.numeric(x) ~ as.factor(cell_type))$p.value})
#pdf('~/elif_ires/FIGURES/KruskalWallis_pvalHistogram_GenesAcrossTissues.pdf', width=5, height=5)
hist(kruskal_pvals,20, xlim = c(0,1))
#dev.off()
length(which(kruskal_pvals < 0.05))
length(kruskal_pvals)
# > 255/278
# [1] 0.9172662
# a1 = aov(as.numeric(dat[5,-1]) ~ as.factor(cell_type))
aov_pvals = apply ( dat[,-1], 1, function (x) {summary.aov(aov(as.numeric(x) ~ as.factor(cell_type)))[[1]][1,5] } )
hist(aov_pvals,20, xlim = c(0,1))
length(which(aov_pvals < 0.01))
length(aov_pvals)
# summary.aov ( a1)
#> length(which(aov_pvals < 0.01))
#[1] 256
#> length(aov_pvals)
#[1] 278
# TukeyHSD ( a1)
# boxplot (as.numeric(dat[5,-1]) ~ as.factor(cell_type))
#pdf ('~/elif_ires/FIGURES/MDS_Replicate_similarity_tissues.pdf', width=5, height=5)
plotMDS(dat[,-1])
#dev.off()

### Define Tissue specificity of IRES
# The idea is to calculate ecdf of each cell type
# Then identify gene entropy; non-uniform distributions will have lower entropy
Mean.IRES = data.frame(ID = dat[,1], ESC.Mean = apply (dat[,2:3], 1, mean), EB.Mean = apply (dat[,4:7], 1, mean), 
            NSC.Mean =  apply (dat[,8:9], 1, mean), Neuron.Mean =  apply (dat[,10:11], 1, mean),
            Limb.Mean =  apply (dat[,12:15], 1, mean)
)

# DENSITIES
par (mfrow = c(2,3))
plot(density(Mean.IRES$ESC.Mean))
plot(density(Mean.IRES$EB.Mean))
plotMDS(Mean.IRES[,-1] )
plot(density(Mean.IRES$NSC.Mean))
plot(density(Mean.IRES$Neuron.Mean))
plot(density(Mean.IRES$Limb.Mean))


my.ecdf = function(x) {ecdf(x)(x)}
Mean.IRES[,2:6] = apply(Mean.IRES[,2:6], 2, my.ecdf)

heatmap.2 (cexCol=.5, as.matrix(Mean.IRES[,-1]), col=redgreen(75), 
           density.info="none", dendrogram="none", 
           scale="none", labRow=F, trace="none" )

Mean.IRES_mean = apply (Mean.IRES[,-1], 1 , mean)
Mean.IRES_centered = Mean.IRES[,-1] - Mean.IRES_mean
Mean.IRES_centered[,6] = Mean.IRES_mean

a1 = apcluster(negDistMat(r=2), as.matrix(Mean.IRES[,-1]))
a2 = apcluster(negDistMat(r=2), Mean.IRES_centered)


######### Test which genes have higher activity than EMCV or HCV in hek_ires dataset
hek_ires = read.csv('~/elif_ires/Elif_DataFiles/072214_HEK_allreplicates_deleted_rows.csv')
hek_ires[,-1] = log10(hek_ires[,-1])

emcv_hek = as.numeric(hek_ires[251,-1])
hcv_hek = as.numeric(hek_ires[97, -1])
prf_hek = as.numeric(hek_ires[250,-1])

test_emcv = function (y){ t.test( as.numeric(y), emcv_hek, alternative="greater")$p.value }
emcv_pvalue = apply ( hek_ires[,-1], 1, test_emcv )
hist ( emcv_pvalue)
length ( which(emcv_pvalue < .05))
length(emcv_pvalue )
# 113 / 279 =  0.4050179

test_hcv = function (y){ t.test( as.numeric(y), hcv_hek, alternative="greater")$p.value }
hcv_pvalue =  apply ( hek_ires[,-1], 1, test_hcv )
hist( hcv_pvalue)
length( which(hcv_pvalue < .05) )
length(hcv_pvalue)
# 75/ 279 = [1] 0.2688172 GENE > HCV

test_prf = function (y){ t.test( as.numeric(y), prf_hek, alternative="greater")$p.value }
prf_pvalue =  apply ( hek_ires[,-1], 1, test_prf )
hist( prf_pvalue)
length( which(prf_pvalue < .05) )
# 275 / 278  ~ .99

stripchart (apply(hek_ires[,-1],1,mean), vertical="T", method= "jitter")


### CALCULATE KAPPA SCORES FOR THE GO TERMS
go_dag = read.table('~/elif_ires/Elif_DataFiles/0728144_Funcassociate_allMGIIDs_tested.txt_Kappa_Network.sif', header=T)
# Calculates kappa similarity between two binary vectors 
calculate_kappa <- function (a1, a2) { 
  Pr_a = sum (!xor(a1,a2)) / length(a1)
  a1_1 = sum(a1) / length(a1)
  a2_1 = sum(a2) /length(a2)
  Pr_e = a1_1 * a2_1 + (1-a1_1) * (1-a2_1)
  kappa = (Pr_a - Pr_e) / (1-Pr_e)
  return (signif(kappa, 2) ) 
}


first_kappas <- c()
for ( i in 1:dim(go_dag)[1]) { 
  for ( j in (i+1):dim(go_dag)[1]) { 
    first_kappas = c(first_kappas, 
                     paste(go_dag[i, 1], go_dag[j, 1],
                           calculate_kappa (go_dag[i, 2:dim(go_dag)[2]], go_dag[j, 2:dim(go_dag)[2]] ), sep="\t"  ) 
    ) 
  }
}

# i = 4 ; j=199 ==> Last calculated j =198
save (first_kappas, file="first_kappas_i4_j198")
load (file="first_kappas_i4_j198")
#write(first_kappas, file = paste(go_dag_joint, "modified", sep="_"))




### OLD PLOS, charts, functions
na.dist <- function(x,...) {
  t.dist <- dist(x,...)
  t.dist <- as.matrix(t.dist)
  t.limit <- 1.1*max(t.dist,na.rm=T)
  t.dist[is.na(t.dist)] <- t.limit
  t.dist <- as.dist(t.dist)
  return(t.dist)
}





