source('~/cancer_jungla/mapCode-colors.r')
library("gplots")


rp_encode_expression = read.table('~/elif/WORM/RP_Expression_Celegans', 
                                  header = T, stringsAsFactors=F, sep="\t")
rp_encode_expression = rp_encode_expression[-grep('rpl-42', rp_encode_expression$Gene.Name), ]
str(rp_encode_expression)
# Remove 'rpl-42 not expressed in most 
mitos = grep('mrp', rp_encode_expression$Gene.Name)
mito_rp_exp = as.matrix(rp_encode_expression[mitos,-c(1:2)])
rp_matrix = as.matrix(rp_encode_expression[-mitos,-c(1:2)])

pdf('~/elif/WORM/FIGURES/RP_Expression_Celegans_ModEncode_log10FPKM.pdf', height=8, width=8)
heatmap.2(log10(rp_matrix), labRow=rp_encode_expression[-mitos,2], scale="none",
          breaks= 50, col = function(x) {colorpanel(x, 'blue3','white', 'red2')},
          cexCol=.75,density.info="none", trace = "none", offsetCol= 0, srtCol=45)
dev.off()

heatmap.2(log10(mito_rp_exp), labRow=rp_encode_expression[mitos,2], scale="none",
          breaks= 50, col = function(x) {colorpanel(x, 'blue3','white', 'red2')},
          cexCol=.75,density.info="none", trace = "none")
