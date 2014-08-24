library(edgeR)
library (gplots)
library(apcluster)
library(hash)
library(xlsx)

firefly = read.xlsx(file='Elif_DataFiles/080814_comparison_raw.xlsx', sheetName="firefly")
renilla = read.xlsx(file='Elif_DataFiles/080814_comparison_raw.xlsx', sheetName="renilla")
colnames(firefly)
colnames(renilla)
firefly[,1] = toupper(firefly[,1])
renilla[,1] = toupper(renilla[,1])

# dat = read.csv ('~/elif_ires/Elif_DataFiles/072314_Elif_comparison_ALLDATA.csv',stringsAsFactors=F)
# dat[,-1] = log10(dat[,-1])
# dat[,1] = toupper(dat[,1])
# Mean.IRES = data.frame(ID = dat[,1], ESC.Mean = apply (dat[,2:3], 1, mean), EB.Mean = apply (dat[,4:7], 1, mean), 
#                        NSC.Mean =  apply (dat[,8:9], 1, mean), Neuron.Mean =  apply (dat[,10:11], 1, mean),
#                        Limb.Mean =  apply (dat[,12:15], 1, mean)
# )

## REmove renilla < 300 from both data sets
# check missing values are overlapping
a2 = which (firefly[,-1] == -1)
a1 = which(renilla[,-1] ==-1)
setdiff(a1,a2)

renilla_threshold = 200 
length ( which ( renilla[,-1] == -1  ) ) / (dim (renilla[,-1])[1] *dim (renilla[,-1])[2])
length ( which ( renilla[,-1] < renilla_threshold ))/ (dim (renilla[,-1])[1] *dim (renilla[,-1])[2])

# Based on http://nar.oxfordjournals.org/content/32/20/e160.full#disp-formula-2
# we define outliers as 
# outlier_detect<- function(r) {
#   range = c(median(r, na.rm=T) + 1.5*IQR(r, na.rm=T), median(r, na.rm=T) - 1.5*IQR(r, na.rm=T))
#   outliers= !(r > range[1] | r < range[2])
#   return(outliers)
# }

firefly[,-1] [renilla[, -1 ] < renilla_threshold] = NA
renilla[,-1] [renilla[, -1 ] < renilla_threshold] = NA

ratios = matrix (nrow = dim (renilla[,-1])[1], ncol = dim (renilla[,-1])[2] )
for ( i in 1: nrow(ratios)) { 
  for ( j in 1:ncol(ratios)){ 
    ratios[i,j]  = firefly[,-1][i,j] / renilla[,-1][i,j]
  }
}
colnames(ratios) = colnames(firefly)[-1]

#write.csv (file ="~/elif_ires/082314_FireflyComparison.csv", firefly,row.names=F )
#write.csv (file ="~/elif_ires/082314_RenillaComparison.csv", renilla,row.names=F )
#write.csv (file = "~/elif_ires/082314_ratios.csv", ratios,row.names=F )

allcors = cor ( ratios[,-c(5,6,14,15)], method = "spearman", use = "pairwise.complete.obs")
dissimilarity <- 1 - cor(allcors)
distance <- as.dist(dissimilarity)
pdf('~/elif_ires/FIGURES/Replicate_hierarchicalCluster_completelinkage_correlation.pdf', height=5, width=5)
plot(hclust(distance, method="complete"))
dev.off()
colSums(is.na(ratios)) / 288

# Decided remove two replicates from Neuron because of clusterin 3-4;
# NSC1-2 > 45% NA so removed

Mean_ratios = data.frame(ID = renilla[,1], 
                         ESC.Mean = log10(apply (ratios[,1:4], 1, mean, na.rm=T) ), 
                         NSC.Mean = log10(apply (ratios[,7:12], 1, mean, na.rm=T) ), 
                         NEU.Mean =  log10(apply (ratios[,c(13,16:18)], 1, mean, na.rm=T) ), 
                         ML.Mean =  log10(apply (ratios[,19:23], 1, mean, na.rm=T) ),
                         EB.Mean =  log10(apply (ratios[,24:27], 1, mean, na.rm=T) )
)
quantile(apply(is.na(Mean_ratios), 1, sum), seq ( 0,1,.1) )
Mean_ratios_complete = Mean_ratios[apply(is.na(Mean_ratios), 1, sum) == 0, ]
dim (Mean_ratios_complete )

par ( las = 1)
par (mfrow = c(2,3))
plot(density(Mean_ratios_complete$ESC.Mean), main = "ESC")
abline (v = Mean_ratios_complete$ESC.Mean[228], col = "red" )
abline ( v  = Mean_ratios_complete$ESC.Mean[2], col = "blue")
plot(density(Mean_ratios_complete$EB.Mean), main = "EB")
abline (v = Mean_ratios_complete$EB.Mean[228], col = "red" )
abline ( v  = Mean_ratios_complete$EB.Mean[2], col = "blue")
plot(density(Mean_ratios_complete$NSC.Mean), main = "NSC")
abline (v = Mean_ratios_complete$NSC.Mean[228], col = "red" )
abline ( v  = Mean_ratios_complete$NSC.Mean[2], col = "blue")
plot(density(Mean_ratios_complete$NEU.Mean), main = "Neuron")
abline (v = Mean_ratios_complete$NEU.Mean[228], col = "red" )
abline ( v  = Mean_ratios_complete$NEU.Mean[2], col = "blue")
plot(density(Mean_ratios_complete$ML.Mean), main = "Mesenchyme")
abline (v = Mean_ratios_complete$ML.Mean[228], col = "red" )
abline ( v  = Mean_ratios_complete$ML.Mean[2], col = "blue")

my.ecdf = function(x) {ecdf(x)(x)}
Mean_ratios_complete[,2:6] = apply(Mean_ratios_complete[,2:6], 2, my.ecdf)

emcv_all = Mean_ratios_complete[228,]
compare_to_emcv = function (x)  {
  emcv_comp = x > emcv_all[-1]
  if (any(emcv_comp)) {
    return (TRUE)
  }
  else {
    return (FALSE)
  }
}

emcv_at_least_one = apply(Mean_ratios_complete[,-1],1,compare_to_emcv)
sum(emcv_at_least_one)

cv_ratios <- apply (Mean_ratios_complete[emcv_at_least_one,-1], 1 , function(x){sd(x)/ mean(x)})
cv.25 = which(cv_ratios > .1)
plot (rowSums(Mean_ratios_complete[emcv_at_least_one,-1])/5 , cv_ratios)
pdf('~/elif_ires/FIGURES/082314_celltype_Gene_CV.1.pdf', width=8, height=8)
h1 = heatmap.2 (cexCol=.5, as.matrix(Mean_ratios_complete[emcv_at_least_one,-1][cv.25,] ), col=redgreen(75), 
                density.info="none", dendrogram="none", 
                scale="none", labRow=Mean_ratios_complete[emcv_at_least_one,][cv.25,1], trace="none", cexRow =.5 )
dev.off()
### CALCULATE MEDIAN IRES FROM THE MEANS
go_dag = readLines('~/elif_ires/Elif_DataFiles/funcassociate_go_associations_mgisymbol.txt')
GO = hash()
GO_full = hash()
# LAST LINE IS EMPTY
# 3910 GO terms has at least one gene
# More than 2; 1680
for (line in go_dag) {
  lineelements = unlist(strsplit(line,split="\t")[[1]])
  GOterm = lineelements[1]
  Genes  = unlist(strsplit(lineelements[3],split=" ")[[1]])
  genes_of_interest = intersect(Genes , Mean_ratios_complete[,1])
  if (length (genes_of_interest) ) {
    GO[[GOterm]] = genes_of_interest
    GO_full[[GOterm]] = lineelements[2]
  }
}

GO_ids = c()
GO_medians = matrix (nrow = length(GO), ncol = 5)
colnames(GO_medians) = c("ESC.Median",  "NSC.Median", "NEU.Median", "ML.Median","EB.Median" )
i = 1
number_of_genes = c()
# Test categories for difference using all genes in addition to taking media
for (key in keys(GO)) {
  if (length(GO[[key]]) > 3 && length(GO[[key]]) < 110) {
    number_of_genes = c(number_of_genes, length(GO[[key]]))  
    GO_medians[i, ] = apply(Mean_ratios_complete[Mean_ratios_complete[,1] %in% GO[[key]],-1] , 2, median)
    GO_ids= c(GO_ids, key)
    i = i+ 1
  }
}
hist(number_of_genes, 50)
quantile(number_of_genes, seq(0,1,.05))
# Cluster GO categories by median
GO_medians = GO_medians[1:length(number_of_genes),]

emcv_atleast_one_logical = apply(GO_medians,1,compare_to_emcv)
sum(emcv_atleast_one_logical)
GO_medians[emcv_atleast_one_logical,]

a2 = apcluster(negDistMat(r=2), GO_medians[emcv_atleast_one_logical,])
GO_medians[emcv_atleast_one_logical,][a2@exemplars,]
for (i in GO_ids[emcv_atleast_one_logical][a2@exemplars]) {
  print(GO_full[[i]])
}

for (i in 1:length(a2@clusters)) {
  print(paste ("Cluster" , i, sep="_" ) )
  for (j in GO_ids[emcv_atleast_one_logical][a2@clusters[[i]] ] ) {
    print (GO_full[[j]] )
  }
  for (k in  a2@clusters[[i]] ) {
    print (GO_medians[emcv_atleast_one_logical,][k,])
  }
}

full_names = c()
for (i in GO_ids[emcv_atleast_one_logical]) {
  full_names = c(full_names, GO_full[[i]])
}

pdf('~/elif_ires/FIGURES/082314_celltype_GO.pdf', width=8, height=8)
h1 = heatmap.2 (cexCol=.5, GO_medians[emcv_atleast_one_logical,], col=redgreen(75), 
                density.info="none", dendrogram="none", 
                scale="none", labRow=full_names, trace="none", cexRow =.2 )
dev.off()

# pdf('~/elif_ires/FIGURES/082314_celltype_GO_10genes.pdf', width=8, height=8)
# h1 = heatmap.2 (cexCol=.5, GO_medians[emcv_atleast_one_logical,], col=redgreen(75), 
#                 density.info="none", dendrogram="none", 
#                 scale="none", labRow=full_names, trace="none", cexRow =.25 )
# dev.off()

cv <- apply (GO_medians[emcv_atleast_one_logical,], 1 , function(x){sd(x)/ mean(x)})
plot(rowSums(GO_medians[emcv_atleast_one_logical,])/5 , cv )
cv_1 = which(cv > .1)

pdf('~/elif_ires/FIGURES/082314_celltype_GO_CV.1.pdf', width=8, height=6)
h1 = heatmap.2 (cexCol=.5, GO_medians[emcv_atleast_one_logical,][cv_1,], col=redgreen(75), 
                density.info="none", dendrogram="row", 
                scale="none", labRow=full_names[cv_1], trace="none", cexRow =.2 )
dev.off()

a3 = apcluster(negDistMat(r=2), GO_medians[emcv_atleast_one_logical,][cv_1,])
GO_medians[emcv_atleast_one_logical,][a3@exemplars,]
for (i in GO_ids[emcv_atleast_one_logical][a3@exemplars]) {
  print(GO_full[[i]])
}

for (i in 1:length(a3@clusters)) {
  print(paste ("Cluster" , i, sep="_" ) )
  for (j in GO_ids[emcv_atleast_one_logical][a3@clusters[[i]] ] ) {
    print (GO_full[[j]] )
  }
  for (k in  a3@clusters[[i]] ) {
    print (GO_medians[emcv_atleast_one_logical,][k,])
  }
}



########## OLD DATA
# test whether mean ratio is the across all cell lines
# We can use either kruskal-wallis non-parametric or aov for parametric assumption
cell_type = c("ESC", "ESC", "EB", "EB","EB","EB", "NSC","NSC", "Neuron","Neuron", "Limb", "Limb" ,"Limb" ,"Limb")
kruskal_pvals = apply ( dat[,-1], 1, function (x) {kruskal.test (as.numeric(x) ~ as.factor(cell_type))$p.value})
#pdf('~/elif_ires/FIGURES/KruskalWallis_pvalHistogram_GenesAcrossTissues.pdf', width=5, height=5)
hist(kruskal_pvals,20, xlim = c(0,1))
#dev.off()
length(which ( p.adjust (kruskal_pvals, method= "fdr" ) < .05 ))
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

# DENSITIES
#pdf ('~/elif_ires/FIGURES/density_plot_replicate_similaritytissues.pdf', width=5, height=5)
par ( las = 2)
par (mfrow = c(2,3))
plot(density(Mean.IRES$ESC.Mean), main = "ESC")
abline (v = log10(2.632078441), col = "red" )
abline ( v  = log10 (3.248532641), col = "blue")
plot(density(Mean.IRES$EB.Mean), main = "EB")
abline (v = log10(0.513328505), col = "red" )
abline ( v  = log10 (0.912015486), col = "blue")
plotMDS(dat[,-1])
plot(density(Mean.IRES$NSC.Mean), main = "NSC")
abline (v = log10(1.223226573), col = "red" )
abline ( v  = log10 (1.094450531), col = "blue")
plot(density(Mean.IRES$Neuron.Mean), main = "Neuron")
abline (v = log10(1.142736044), col = "red" )
abline ( v  = log10 (2.827249104), col = "blue")
plot(density(Mean.IRES$Limb.Mean), main = "Mesenchyme")
abline (v = log10(3.38589839), col = "red" )
abline ( v  = log10 (4.07622128), col = "blue")
#dev.off()


my.ecdf = function(x) {ecdf(x)(x)}
Mean.IRES[,2:6] = apply(Mean.IRES[,2:6], 2, my.ecdf)

pdf ('~/elif_ires/FIGURES/heatmap_tissues_rankscore_normalized.pdf', width=5, height=5)
heatmap.2 (cexCol=.5, as.matrix(Mean.IRES[,-1]), col=redgreen(75), 
           density.info="none", dendrogram="none", 
           scale="none", labRow=F, trace="none" )
dev.off()
# Mean.IRES_mean = apply (Mean.IRES[,-1], 1 , mean)
# Mean.IRES_centered = Mean.IRES[,-1] - Mean.IRES_mean
# Mean.IRES_centered[,6] = Mean.IRES_mean


a1 = apcluster(negDistMat(r=2), as.matrix(Mean.IRES[,-1]))
Mean.IRES[a1@exemplars,]
pdf ('~/elif_ires/FIGURES/affinity_propogationclustering_heatmap.pdf', width=5, height=5)
heatmap(a1)
dev.off()
#plot(a1 , as.matrix(Mean.IRES[,-1]))
# a2 = apcluster(negDistMat(r=2), Mean.IRES_centered)

for ( i in 1: length ( a1@clusters))  { 
  write.csv ( Mean.IRES[a1@clusters[[i]], ] , row.names = F, 
                file = paste ("~/elif_ires/APCLUSTERS/Genes_inCluster_", i, ".xls", sep = "" ) )
}



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





