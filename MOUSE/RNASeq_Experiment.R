read_count_million = 400
luciferase_rpkm = 15
number_of_constructs = 120
luciferase_length = 3000
junc_pos1 = 1000
junc_pos2 = 2000
len_on_each_side = 175

number_of_reads_overlapping_luciferase = (read_count_million * luciferase_rpkm * luciferase_length) / 1000
per_construct = number_of_reads_overlapping_luciferase / number_of_constructs

dist_num1 = dist_num2 = c()
for ( dummy in 1:1000) {
read_starts = sample (1:luciferase_length,size=per_construct,replace=T)  
num1 = 0
num2 = 0
for (i in 1:length(read_starts)) {
  if ( read_starts[i] > junc_pos1 - len_on_each_side  && read_starts[i] < junc_pos1 + len_on_each_side ) {
    num1 = num1 + 1
  }
  if ( read_starts[i] > junc_pos2 - len_on_each_side  && read_starts[i] < junc_pos2 + len_on_each_side ) {
    num2 = num2 + 1
  }
}
dist_num1 = c(dist_num1, num1)
dist_num2 = c(dist_num2, num2)
}
hist(dist_num1)
hist(dist_num2)
hist ( dist_num1 - dist_num2)
quantile ( dist_num1 - dist_num2, seq(0,1,.05))

ratios = c()
for ( dum in 1:1000) { 
different_ratio = sample (c(0,1), prob=c(.4,.6), size = 1000, replace = T)
ratios = c(ratios, sum ( different_ratio) / 1000 ) 
}

remaining = 150 - (per_construct * ratios) 
hist((per_construct * ratios)  - remaining )
