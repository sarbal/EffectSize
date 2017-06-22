load("data/all_gene_sets.Rdata")
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


