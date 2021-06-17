library(ExactMultinom)
library(dplyr, warn.conflicts = F)

args <- commandArgs(trailingOnly = T)
data <- read.csv(file = args[1], header = F, skip = 1, sep=" ")
# print(data)
# apply(data, 1, function(row) { multinom.test(row, rep(1, length(row)), method="Monte-Carlo")["pvals_mc"][[1]]
# })

new_data <- data %>% rowwise() %>% mutate(p = unlist(multinom.test(c_across(), rep(1, ncol(data)), method="Monte-Carlo")["pvals_mc"])[[2]])
print(new_data$p)
