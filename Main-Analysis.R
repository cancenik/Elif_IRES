library(edgeR)
dat = read.csv ('./072314_Elif_comparison_ALLDATA.csv')
dat[,-1] = log10(dat[,-1])

# test whether mean 
cell_type = c("ESC", "ESC", "EB", "EB","EB","EB", "NSC","NSC", "Neuron","Neuron", "Limb", "Limb" ,"Limb" ,"Limb")
kruskal_pvals = apply ( dat[,-1], 1, function (x) {kruskal.test (as.numeric(x) ~ as.factor(cell_type))$p.value})
hist(kruskal_pvals,20, xlim = c(0,1))
length(which(kruskal_pvals < 0.05))
length(kruskal_pvals)
# > 255/278
# [1] 0.9172662
# a1 = aov(as.numeric(dat[5,-1]) ~ as.factor(cell_type))
aov_pvals = apply ( dat[,-1], 1, function (x) {summary.aov(aov(as.numeric(x) ~ as.factor(cell_type)))[[1]][1,5] } )
hist(aov_pvals,20, xlim = c(0,1))
length(which(aov_pvals < 0.05))
length(aov_pvals)
# summary.aov ( a1)
# TukeyHSD ( a1)
# boxplot (as.numeric(dat[5,-1]) ~ as.factor(cell_type))
plotMDS(dat[,-1])
