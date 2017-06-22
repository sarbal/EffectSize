load("data/asd_gene_sets.Rdata")
load("data/all.rank.properties.Rdata")
outdir = "out/"

# Running gene and gene score enrichment tests
teste =calculate_functional_effects_enrichments(asd.genes.tally, func.props.list[[2]], studies, nsub)                       # Gene properties
testd =calculate_functional_effects_enrichments(asd.genes.tally, func.props.list[[4]], studies, nsub)                       # DE


testg = calculate_functional_effects_enrichments_binary(asd.genes.tally, func.props.list[[3]], studies=studies, nsub=nsub)  # Gene sets
testk = calculate_functional_effects_enrichments_binary(asd.genes.tally, func.props.list[[6]],  studies=studies, nsub=nsub) # KEGG
testgo = calculate_functional_effects_enrichments_binary(asd.genes.tally, func.props.list[[5]], studies=studies, nsub=nsub) # GO

# Running network connectivity tests
testn = calculate_functional_effects_networks(asd.genes.tally, net, "network", outdir, studies, nsub)                       # Network connectivity

pvals.all.pre  = cbind( teste[[2]][,1], testd[[2]][,1], as.numeric(testgo[[1]][,5]), as.numeric(testk[[1]][,5]), as.numeric(testg[[1]][,5]), t(testn[[1]])[,3] )
pvals.all.post = cbind( teste[[2]][,2], testd[[2]][,2], as.numeric(testgo[[2]][,5]), as.numeric(testk[[2]][,5]), as.numeric(testg[[2]][,5]), t(testn[[2]])[,3] )

func.labels = c( colnames(func.props.list[[2]]), colnames(func.props.list[[4]]), colnames(func.props.list[[3]]), colnames(func.props.list[[6]]), colnames(func.props.list[[5]]))

# Correlations
slopes = calc_slopes(func.labels, effect.sizes[filt.sub], pvals.all.post[,filt.sub], rr = 1000)
