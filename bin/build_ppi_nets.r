# example
library(EGAD)
load(file="intact.human.Rdata")
data.mod = intact.mod
genes = unique(c(data.mod[,2], data.mod[,1]))
net = build_binary_network(data.mod, genes)
ext_net = extend_network(net)
