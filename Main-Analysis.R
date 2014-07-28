library(edgeR)
dat = read.csv ('~/elif_ires/Elif_DataFiles/072314_Elif_comparison_ALLDATA.csv')
dat[,-1] = log10(dat[,-1])
hek_ires = read.csv('~/elif_ires/Elif_DataFiles/072214_HEK_allreplicates_deleted_rows.csv')
hek_ires[,-1] = log10(hek_ires[,-1])

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


# Test which genes have higher activity than EMCV or HCV in hek_ires dataset
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


stripchart (log10(ratio$ratio), vertical="T", method= "jitter", add = TRUE, at = 0.7))

ires = read.csv('~/Google Drive/barna lab/Experiments/5_UTR_Conservation/Transfection_bicistronic_luciferase /SGI_screen/comparison/042014_comparison_different_cellTypes.csv')
ires_hek = ires[, 1:2]
ires_esc = ires[, 3:4]
ires_limbs = ires[,5:6]
ires_nsc = ires[,7:8]

ires_comparison1 = merge (ires_hek,ires_esc, by.x="DNA_HEK",by.y= "DNA_mESC")

ires_comparison2 = merge (ires_comparison1,ires_limbs, by.x="DNA_HEK",by.y= "DNA_limbs")
ires_comparison3 = merge (ires_comparison2,ires_nsc, by.x="DNA_HEK",by.y= "DNA_NSC")
ires_comparison_final = unique (ires_comparison3, incomparables = FALSE)
morethanoneoccurence <- duplicated(ires_comparison_final[,1])
dim(ires_comparison_final[!morethanoneoccurence, ])
first_selected_ires_comparison = ires_comparison_final[!morethanoneoccurence, ]
first_selected_ires_comparison [,2]  <- as.numeric (as.character (first_selected_ires_comparison [,2]))
library (gplots)
ires.matrix = as.matrix (first_selected_ires_comparison [,2:5])
#ires.matrix[is.na(ires.matrix)] <- -1
na.dist <- function(x,...) {
  t.dist <- dist(x,...)
  t.dist <- as.matrix(t.dist)
  t.limit <- 1.1*max(t.dist,na.rm=T)
  t.dist[is.na(t.dist)] <- t.limit
  t.dist <- as.dist(t.dist)
  return(t.dist)
}
heatmap.2 (ires.matrix, distfun=na.dist, col=redgreen(75), 
           density.info="none", dendrogram="row", 
           scale="row", labRow=F, 
           trace="none", na.color="blue")

nas = apply ( is.na (ires.matrix),1, any )
ires.matrix2 = ires.matrix[!nas, ]
heatmap.2 (ires.matrix2, distfun=na.dist, col=redgreen(75), 
           density.info="none", dendrogram="none", 
           labRow = first_selected_ires_comparison[!nas,1],cexRow=.2,
           scale="row",
           trace="none", na.color="blue")

cor (ires.matrix2, method = "spearman")


### CALCULATE KAPPA SCORES FOR THE GO TERMS
go_dag = go_dag_abs8
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
#write(first_kappas, file = paste(go_dag_joint, "modified", sep="_"))


