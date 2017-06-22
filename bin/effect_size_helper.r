#########################################
#  	Helper functions	    	#
#########################################

##  Written: Sara Ballouz
##  Date: April 28th, 2016

library(MASS)
library(Matrix)
library(scales)
library(gplots)
library(zoo)
library(plyr)
library(RColorBrewer)


# helper functions
residuals <- function(x,y,A,B,C){ (A*x + B*y + C) }

# Transparent colors
makeTransparent<-function(someColor, alpha=100)
{
	newColor<-col2rgb(someColor)
	apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
	blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}



geo_mean <- function(data) {
    log_data <- log(data)
    gm <- exp(mean(log_data[is.finite(log_data)]))
    return(gm)
}


make_net <- function(gene.corr,n){
	#tic()
 print(".")
	net <- matrix(rank(gene.corr, na.last = "keep", ties.method = "average"), nrow = n, ncol = n)
	print(".")
	rownames(net) <- rownames(gene.corr)
	colnames(net) <- colnames(gene.corr)
	net <- net/max(net, na.rm = T)
	print(".")
	diag(net) <- 1
	print(".")
	#toc()
	return(net)
}


# sub sampling functions
sub_sample <- function(X, n, r){
	if( length(X) < n) {
		mean(X, na.rm=T)
	} else {
		mean(sapply( 1:r, function(i) mean(sample( X, n, replace=F),na.rm=T) ))

	}
}


sub_sample_cond <- function(X,n, r, g.l.c){
	if( length(X) < n) {
		mean(X, na.rm=T)
	} else {
                sample = sample( unique(g.l.c[,2]) ,n) # pick locus/condition
		nsub.samp = sapply(1:n, function(i) sample( rep(which( g.l.c[,2] == sample[i]),2) ,1)) #pick "gene" from locus/condition
		mean(sapply( 1:r, function(i) mean( X[nsub.samp] ,na.rm=T) ))
	}
}


sub_sample_enrich <- function(X, g.l, v, n, r ){
        print(length(X))
	if( length(X) <= n) {
		gene_set_enrichment(X, g.l, v)
	} else {
		gene_set_enrichment_sample2( X,g.l, v , n, r )

	}
}



sub_sample_enrich_cond <- function(X, g.l, v, n, r, g.l.c ){
        print(length(X))
	if( length(X) <= n) {
		gene_set_enrichment(X, g.l, v)
	} else {
		gene_set_enrichment_sample_cond( X,g.l, g.l.c, v , n, r )

	}
}



sub_sample_resd <- function(net, dir, label, g.l, n, r, f ){

        if( dim(g.l)[1] < n) {
                #c(0,0,1)
                residual_connectivity_score(net,dir,label,g.l)
	} else {
		residual_connectivity_score_resampling(net, dir, label, g.l, n, r, f)

	}
}



sub_sample_resd_cond <- function(net, dir, label, g.l, n, r, f, g.c ){

        if( dim(g.l)[1] < n) {
                #c(0,0,1)
                residual_connectivity_score(net,dir,label,g.l)
	} else {
		residual_connectivity_score_resampling_conditional(net, dir, label, g.l, n, r, f, g.c)

	}
}


sub_sample_wilcox <- function(X.sub, X.all, n, r, alt ){
	if(length(X.sub) <=0){
		return(1)
	}
	else if( length(X.sub) <= n) {
		wilcox.test(X.sub, X.all, alt=alt)$p.value
	} else {
                geo_mean(sapply(1:r, function(i) wilcox.test(sample(X.sub, n, replace=F), X.all, alt=alt)$p.value))
	}
}


sub_sample_wilcox_cond <- function(X.sub, X.all, n, r, alt, g.l.c ){
	if(length(X.sub) <=0){
		return(1)
	}
	else if( length(X.sub) <= n) {
		wilcox.test(X.sub, X.all, alt=alt)$p.value
	} else {
                sample = sample( unique(g.l.c[,2]) ,n) # pick locus/condition
		nsub.samp = sapply(1:n, function(i) sample( rep(which( g.l.c[,2] == sample[i]),2) ,1)) #pick "gene" from locus/condition

                geo_mean(sapply(1:r, function(i) wilcox.test(X.sub[nsub.samp], X.all, alt=alt)$p.value))
	}
}


# Property tests
 


### Gene set enrichment tests
gene_set_enrichment <- function(genes, genes.labels, voc){

        genes.names = rownames(genes.labels)
	labels.names = colnames(genes.labels)
        genes.counts = rowSums(genes.labels)
	labels.counts = colSums(genes.labels)              			# p

	m = match ( genes, genes.names )
	filt.genes  = !is.na(m)
	filt.labels = m[filt.genes]


	labels.counts.set = rep( sum(filt.genes), length(labels.counts) )	# g

        m = match (labels.names, voc[,1])
        v.f = !is.na(m)
        v.g = m[v.f]

	universe = rep ( dim(genes.labels)[1], dim(genes.labels)[2])
	if(  length(filt.labels) == 1 ) { genes.counts.set = genes.labels[filt.labels,] }
	else { genes.counts.set = colSums(genes.labels[filt.labels,]) }             ## does weird things with 0 sets

	test =  cbind( (genes.counts.set -1) , labels.counts, universe-labels.counts, labels.counts.set)
        pvals = phyper(test[,1], test[,2], test[,3], test[,4], lower.tail=F)
        sigs = pvals < ( 0.05/length(pvals) )
        pvals.adj = p.adjust( pvals, method="BH")

	results = cbind(voc[v.g,1:2], test[v.f,c(1,2)], pvals[v.f], pvals.adj[v.f], sigs[v.f])
	return (results)

}


gene_set_enrichment_sample <- function(genes, genes.labels, voc, n, r){

        genes.names = rownames(genes.labels)
	labels.names = colnames(genes.labels)
        genes.counts = rowSums(genes.labels)
	labels.counts = colSums(genes.labels)              			# p
        m = match (labels.names, voc[,1])
        v.f = !is.na(m)
        v.g = m[v.f]
        universe = rep ( dim(genes.labels)[1], dim(genes.labels)[2])
        pvals.sum = rep(0, sum(v.f))
        pvals.adj.sum = pvals.sum
        for (i in 1:r){
                sub = sample(genes, n, replace=F)
                m = match ( sub, genes.names )
                filt.genes  = !is.na(m)
                filt.labels = m[filt.genes]
                labels.counts.set = rep( sum(filt.genes), length(labels.counts) )	# g

                if(  length(filt.labels) == 1 ) { genes.counts.set = genes.labels[filt.labels,] }
                else { genes.counts.set = colSums(genes.labels[filt.labels,]) }             ## does weird things with 0 sets

                test =  cbind( (genes.counts.set -1) , labels.counts, universe-labels.counts, labels.counts.set)
                pvals = phyper(test[,1], test[,2], test[,3], test[,4], lower.tail=F)
                pvals.adj = p.adjust( pvals, method="BH")
                pvals.sum = pvals  + pvals.sum
                pvals.adj.sum = pvals.adj + pvals.adj.sum
        }
        pvals = pvals.sum/r
        pvals.adj = pvals.adj.sum/r
        sigs = pvals < ( 0.05/length(pvals) )


        results = cbind(voc[v.g,1:2], test[v.f,c(1,2)], pvals[v.f], pvals.adj[v.f], sigs[v.f])
	return (results)

}



gene_set_enrichment_sample2 <- function(genes, genes.labels, voc, n, r){

        genes.names = rownames(genes.labels)
	labels.names = colnames(genes.labels)
        genes.counts = rowSums(genes.labels)
	labels.counts = colSums(genes.labels)              			# p
        m = match (labels.names, voc[,1])
        v.f = !is.na(m)
        v.g = m[v.f]
        universe = rep ( dim(genes.labels)[1], dim(genes.labels)[2])
        pvals.sum = matrix(0, nrow=length(labels.names), ncol=r)
        pvals.adj.sum = pvals.sum
        for (i in 1:r){
                sub = sample(genes, n, replace=F)
                m = match ( sub, genes.names )
                filt.genes  = !is.na(m)
                filt.labels = m[filt.genes]
                labels.counts.set = rep( sum(filt.genes), length(labels.counts) )	# g

                if(  length(filt.labels) == 1 ) { genes.counts.set = genes.labels[filt.labels,] }
                else { genes.counts.set = colSums(genes.labels[filt.labels,]) }             ## does weird things with 0 sets

                test =  cbind( (genes.counts.set -1) , labels.counts, universe-labels.counts, labels.counts.set)
                pvals = phyper(test[,1], test[,2], test[,3], test[,4], lower.tail=F)
                pvals.adj = p.adjust( pvals, method="BH")
                pvals.sum[,i] = pvals
                pvals.adj.sum[,i] = pvals.adj
        }
        pvals = apply( pvals.sum, 1, geo_mean)
        pvals.adj = apply( pvals.adj.sum, 1, geo_mean)
        sigs = pvals.adj < ( 0.05/length(pvals.adj) )


        results = cbind(voc[v.g,1:2], test[v.f,c(1,2)], pvals[v.f], pvals.adj[v.f], sigs[v.f])
	return (results)

}



#### Network connectivity tests

residual_connectivity_score <- function(network, dir, label, gene.list){
        if( length(gene.list) < 0 ){ return (list(0,0,1)) }
	m = match ( rownames(network), attr$entrezID)
	f = !is.na(m)
	g = m[f]

	m = match( gene.list[,1], attr$name[g] )
	f.r = !is.na(m)
	f.ar = m[f.r]

        if( sum(f.r) <= 2 ) { return (list(0,0,1)) }

	network_sub = network[f,f][f.ar,f.ar]

	node_degree =  rowSums(network, na.rm=T)
        print( mean(node_degree, na.rm=T) )
	print ( sd(node_degree, na.rm=T) )

	node_degree_sub =  rowSums(network_sub, na.rm=T)
	print( mean(node_degree_sub, na.rm=T) )
	print ( sd(node_degree_sub, na.rm=T) )

	N = dim(network)[1]
	n = dim(network_sub)[1]
	x = c(0,N)
	y = c( n/N * x[1], n )
	print(N)
	print(n)
        if( n <= 2 ) { return (list(0,0,1)) }
	#h.all = get_density(hist(node_degree, plot=F))
	#h.sub = get_density(hist(node_degree[f][f.ar], plot=F))

	resd.sub = residuals( node_degree[f][f.ar], node_degree_sub, -length( node_degree_sub), N, 0 )

#	points = gene.list[f.r,2]
#	p      = gene.list[f.r,2] > 1
#	genesp = gene.list[f.r,1]

	resd.rand  = unlist(explore_sub_network_list_random3(network, dir, paste("random",label, sep="."),n))
	resd.rand = as.matrix(resd.rand )
	resd.rand = sort(resd.rand)
#	h.rand = get_density(hist(resd.rand , plot=F))

	X = cbind(node_degree[f][f.ar], node_degree_sub)
	H = X %*% solve( t(X) %*% X ) %*% t(X)
	X = resd.sub/ sd(resd.sub) * sqrt(1 - diag(H) )

     	a=round(mean(resd.rand , na.rm=T),3)
#	b=round(mean(X[p] , na.rm=T),3)
#	c=round(mean(X[!p] , na.rm=T),3)
        d=round(mean(X , na.rm=T),3)
        test = wilcox.test(resd.sub, resd.rand, alt="g")

	pvals = sapply( 1:length(X), function(i) sum ( resd.rand > X[i] )  )  / length(resd.rand)
	pvals.adj = p.adjust( pvals, method="BH")
	print(c(a, d, test$p.value))

        results.all = cbind( as.character(attr$name[g][f.ar]),  as.character(attr$entrezID[g][f.ar]), node_degree[f][f.ar], node_degree_sub, resd.sub, X, pvals, pvals.adj)
	colnames(results.all) = c("Gene", "entrezID", "Node degree full", "Node degree sub", "Residuals", "X", "P-vals", "P-vals adj")
        write.table(results.all, file = paste(dir, label, ".residuals", sep=""))
	write.table(c(N,n), file = paste(dir, label, ".Ns", sep=""))

        #return(results.all)
        return (list(a, d, test$p.value))

}



residual_connectivity_score_resampling <- function(network, dir, label, gene.list, nsub, r, file){
        if( length(gene.list) < 0 ){ return (list(0,0,1)) }
	m = match ( rownames(network), attr$entrezID)
	f = !is.na(m)
	g = m[f]

	m = match( gene.list[,1], attr$name[g] )
	f.r = !is.na(m)
	f.ar = m[f.r]

        if( sum(f.r) <= nsub ) { return (list(0,0,1)) }

        network_sub = network[f,f][f.ar,f.ar]

	node_degree =  rowSums(network, na.rm=T)
        print( mean(node_degree, na.rm=T) )
	print ( sd(node_degree, na.rm=T) )

	node_degree_sub =  rowSums(network_sub, na.rm=T)
	print( mean(node_degree_sub, na.rm=T) )
	print ( sd(node_degree_sub, na.rm=T) )

	N = dim(network)[1]
	n = dim(network_sub)[1]
	x = c(0,N)
	y = c( n/N * x[1], n )
	print(N)
	print(n)
        if( n <= nsub ) { return (list(0,0,1)) }
	#h.all = get_density(hist(node_degree, plot=F))
	#h.sub = get_density(hist(node_degree[f][f.ar], plot=F))


        if( file.exists( file )){
                resd.rand = unlist(read.table(file))
	} else {
                resd.rand  = unlist(explore_sub_network_list_random3(network, dir, paste("random.sub.",label, sep="."),nsub))
        }

        resd.rand = as.matrix(resd.rand )
	resd.rand = sort(resd.rand)

        tic();
        mat=matrix(0, ncol=5, nrow=r)
        for( k in 1:r){
                nsub.samp = sample(n, nsub, replace=F)
        	node_degree_sub =  rowSums(network_sub[nsub.samp,nsub.samp])

                resd.sub = residuals( node_degree[f][f.ar][nsub.samp], node_degree_sub, -length( node_degree_sub), N, 0 )

		# studentized residuals
                X = cbind(node_degree[f][f.ar][nsub.samp], node_degree_sub)
                H = X %*% solve( t(X) %*% X ) %*% t(X)
                X = resd.sub/ sd(resd.sub) * sqrt(1 - diag(H) )

                pvals = sapply( 1:length(X), function(i) sum ( resd.rand > X[i] )  )  / length(resd.rand)
                pvals.adj = p.adjust( pvals, method="BH")
                a = round(mean(resd.sub , na.rm=T),3)
                b = mean( pvals, na.rm=T)
                c = mean( pvals.adj, na.rm=T)
                d=round(mean(X , na.rm=T),3)
                test = wilcox.test(resd.sub, resd.rand, alt="g")
                mat[k,] = c(a,b,c,d, test$p.value)
                # print(c(a,b,c,d, test$p.value))

         }
         toc();
          write.table(mat, file = paste(dir, label, ".resamp.residuals", sep=""))

         a = round(mean(resd.rand , na.rm=T),3)
         d = round( mean(mat[,4], na.rm=T),3)

         test$p.value = geo_mean(mat[,5])
         # test$p.value = mean(mat[,5], na.rm=T)

#        results.all = cbind( as.character(attr$name[g][f.ar]),  as.character(attr$entrezID[g][f.ar]), node_degree[f][f.ar], node_degree_sub, resd.sub, X, pvals, pvals.adj)
#	colnames(results.all) = c("Gene", "entrezID", "Node degree full", "Node degree sub", "Residuals", "X", "P-vals", "P-vals adj")
#        write.table(results.all, file = paste(dir, label, ".residuals", sep=""))
	#return(results.all)
        return (list(a, d, test$p.value))

}



residual_connectivity_score_resampling_conditional <- function(network, dir, label, gene.list, nsub, r, file, gene.list.conditions){
        if( length(gene.list) < 0 ){ return (list(0,0,1)) }
	m = match ( rownames(network), attr$entrezID)
	f = !is.na(m)
	g = m[f]

	m = match( gene.list[,1], attr$name[g] )
	f.r = !is.na(m)
	f.ar = m[f.r]

        if( sum(f.r) <= nsub ) { return (list(0,0,1)) }

        network_sub = network[f,f][f.ar,f.ar]

	node_degree =  rowSums(network, na.rm=T)
        print( mean(node_degree, na.rm=T) )
	print ( sd(node_degree, na.rm=T) )

	node_degree_sub =  rowSums(network_sub, na.rm=T)
	print( mean(node_degree_sub, na.rm=T) )
	print ( sd(node_degree_sub, na.rm=T) )

	N = dim(network)[1]
	n = dim(network_sub)[1]
	x = c(0,N)
	y = c( n/N * x[1], n )
	print(N)
	print(n)
        if( n <= nsub ) { return (list(0,0,1)) }
	#h.all = get_density(hist(node_degree, plot=F))
	#h.sub = get_density(hist(node_degree[f][f.ar], plot=F))


        if( file.exists( file )){
                resd.rand = unlist(read.table(file))
	} else {
                resd.rand  = unlist(explore_sub_network_list_random3(network, dir, paste("random.sub.",label, sep="."),nsub))
        }

        resd.rand = as.matrix(resd.rand )
	resd.rand = sort(resd.rand)

        tic();
        mat=matrix(0, ncol=5, nrow=r)
        for( k in 1:r){
                # nsub.samp = sample(n, nsub, replace=F)
                sample = sample( unique(gene.list.conditions[f.r,2]) ,nsub) # pick locus/condition
                nsub.samp = sapply(1:nsub, function(i) sample( rep(which( gene.list.conditions[f.r,2] == sample[i]),2) ,1)) #pick "gene" from locus/condition

        	node_degree_sub =  rowSums(network_sub[nsub.samp,nsub.samp], na.rm=T)

                resd.sub = residuals( node_degree[f][f.ar][nsub.samp], node_degree_sub, -length( node_degree_sub), N, 0 )

		# studentized residuals
                X = cbind(node_degree[f][f.ar][nsub.samp], node_degree_sub)
                H = X %*% solve( t(X) %*% X ) %*% t(X)
                X = resd.sub/ sd(resd.sub) * sqrt(1 - diag(H) )

                pvals = sapply( 1:length(X), function(i) sum ( resd.rand > X[i] )  )  / length(resd.rand)
                pvals.adj = p.adjust( pvals, method="BH")
                a = round(mean(resd.sub , na.rm=T),3)
                b = mean( pvals, na.rm=T)
                c = mean( pvals.adj, na.rm=T)
                d=round(mean(X , na.rm=T),3)
                test = wilcox.test(resd.sub, resd.rand, alt="g")
                mat[k,] = c(a,b,c,d, test$p.value)
                # print(c(a,b,c,d, test$p.value))

         }
         toc();
          write.table(mat, file = paste(dir, label, ".resamp.residuals", sep=""))

         a = round(mean(resd.rand , na.rm=T),3)
         d = round( mean(mat[,4], na.rm=T),3)

         test$p.value = geo_mean(mat[,5])
         # test$p.value = mean(mat[,5], na.rm=T)

#        results.all = cbind( as.character(attr$name[g][f.ar]),  as.character(attr$entrezID[g][f.ar]), node_degree[f][f.ar], node_degree_sub, resd.sub, X, pvals, pvals.adj)
#	colnames(results.all) = c("Gene", "entrezID", "Node degree full", "Node degree sub", "Residuals", "X", "P-vals", "P-vals adj")
#        write.table(results.all, file = paste(dir, label, ".residuals", sep=""))
	#return(results.all)
        return (list(a, d, test$p.value))

}



explore_sub_network_list_random3 <- function(network, dir, label, nR){
	m = match ( rownames(network), attr$entrezID)
	f = !is.na(m)
	g = m[f]

        node_degree =  rowSums(network, na.rm=T)
        print( mean(node_degree, na.rm=T) )
	print ( sd(node_degree, na.rm=T) )
        filts = !is.na(attr$entrezID[g])
	N = dim(network)[1]
        network = network[f,f][filts,filts]
        summary = matrix(0, ncol=2, nrow=1000)
        data  = list()
        for ( i in 1:1000){

                gene.list = sample(attr$name[g][filts], nR)
        	m = match( gene.list, attr$name[g][filts] )
                f.r = !is.na(m)
                f.ar = m[f.r]

                network_sub = network[f.ar,f.ar]
                node_degree_sub =  rowSums(network_sub, na.rm=T)

                n = dim(network_sub)[1]
        	x = c(0,N)
                y = c( n/N * x[1], n )
		k = 1
		data.sub = list()
		#for( j in c(1,2,5,10) ){
		#	frac = 10/j
			frac = 1
	                prb   = sample( gene.list[f.r], length(gene.list[f.r])/frac )
	                m = match( gene.list[f.r], prb)
	                prb = !is.na(m)
	                network_sub.prb     = network_sub[prb,prb]
	                node_degree_sub.prb = rowSums(network_sub.prb, na.rm=T)
	                X = residuals( node_degree[f][filts][f.ar][prb], node_degree_sub.prb, -length( node_degree_sub.prb), N, 0 )
	                H = X %*% solve( t(X) %*% X ) %*% t(X)
	                resd.sub.prb = X/(sd(X)*sqrt(1-diag(H)))
	                summary[i,k] = mean(resd.sub.prb, na.rm=T)
                	summary[i,k+1] = sd(resd.sub.prb, na.rm=T)
                #	data.sub[[j]] =  resd.sub.prb
                #	k = k + 2
                #}
                data[[i]] = resd.sub.prb
        }
	write.table(summary, file=paste(dir,label,"summary.stats.random", sep="" ) )
	write.table(data, file=paste(dir,label,"residuals.random", sep="" ) )
	return(data)
}





calculate_functional_effects_networks <- function( gene.sets, network, nettype, dir, studies, nsub, r=1000, j =2 ){
        studies.genes = rownames(gene.sets)

	genes = lapply(1:length(studies), function(i) as.matrix( cbind(studies.genes, as.numeric(gene.sets[,i]))[gene.sets[,i]>0,] ) )
	labels =  lapply(1:length(studies), function(i) paste(studies[i] ,nettype,sep=".") )

        file = paste(dir,"random.", nsub, ".",nettype,".residuals.random",sep="" )

        an = lapply(1:length(studies), function(i) residual_connectivity_score(network , dir, labels[[i]], genes[[i]] ) )
	an = matrix(unlist(an), ncol =length(studies), nrow=3, byrow=F)
        dn = lapply(1:length(studies), function(i) sub_sample_resd(network , dir,labels[[i]],genes[[i]], nsub, r, file) )
	dn = matrix(unlist(dn), ncol =length(studies), nrow=3, byrow=F)

	colnames(an) = studies
	rownames(an) = c("res","rand", "p")
	colnames(dn) = studies
	rownames(dn) = c("res","rand", "p")

       return( list( an,dn))
}




calculate_functional_effects_enrichments_binary <- function( gene.sets, func.props, voc, studies, nsub, alt="g", r1=1000, r2=100 ){
        if( missing(voc)) {
		voc = matrix(rep( colnames(func.props), 3 ), byrow=F, ncol=3)
	}

	studies.genes = rownames(gene.sets)
        props.genes = rownames(func.props)
	res.test = list()
	pvals.test = list()

	res     = lapply(1:length(studies), function(i) gene_set_enrichment( studies.genes[gene.sets[,i]>0], func.props, voc ))
        res.sub = lapply(1:length(studies), function(i) sub_sample_enrich( studies.genes[gene.sets[,i]>0], func.props, voc, nsub,r1) )

        pvals     = sapply(1:length(studies), function(i)  res[[i]][,6] )
        pvals.sub = sapply(1:length(studies), function(i)  res.sub[[i]][,6] )

      return(list(res, res.sub, pvals, pvals.sub) )
}


calculate_functional_effects_enrichments <- function( gene.sets, func.props, studies, nsub, alt="g", r1=1000, r2=100  ){
        studies.genes = rownames(gene.sets)
        props.genes = rownames(func.props)
 	res.test = list()
	pvals.test = list()
	for( I in 1:dim(func.props)[2]){
		f.n = !is.na(func.props[,I])
		m = match( studies.genes, props.genes[f.n])
		f.ag = !is.na(m)
		f.g = m[f.ag]

                prop = as.numeric(func.props[f.n,I])

                res     = sapply(1:length(studies), function(i) mean( prop[f.g][gene.sets[f.ag,i]>0], na.rm=T))
                res.sub = sapply(1:length(studies), function(i) sub_sample( prop[f.g][gene.sets[f.ag,i]>0],nsub,r1) )

                pvals     = sapply(1:length(studies), function(i) if( length(prop[f.g][gene.sets[f.ag,i]>0]) > 0) { wilcox.test(prop[f.g][gene.sets[f.ag,i]>0], prop , alt=alt)$p.value} else { return(1)}  )
                pvals.sub = sapply(1:length(studies), function(i) sub_sample_wilcox( prop[f.g][gene.sets[f.ag,i]>0],prop,nsub,r2, alt))

		res.test[[I]] = cbind(res,res.sub)
		pvals.test[[I]] = cbind(pvals, pvals.sub)
	}

       return(list(res.test, pvals.test) )

}




calc_slopes <- function(func.labels,effect.sizes, pvals.all, PERMUTE_FLAG=TRUE, rr=100){
	labels = c("Raw(pearson corr)", "Ranked(spearman corr)", "Ranked(ties averaged)", "Ranked(ordered)", "Slopes(ordered)", "Slopes(ranked)") #, "Raw (slopes)")

	rank.order = 1:length(effect.sizes)
	rank.effect = rank(effect.sizes)
	cor.p = sapply(1:length(func.labels), function(i) cor( effect.sizes, -log10(pvals.all[,i] ), method="p") )
	cor.s = sapply(1:length(func.labels), function(i) cor( effect.sizes, -log10(pvals.all[,i] ), method="s") )

	cor.rs = sapply(1:length(func.labels), function(i) cor( rank.order, rank(-log10(pvals.all[,i]))))
	cor.r  = sapply(1:length(func.labels), function(i) cor( rank.effect, rank(-log10(pvals.all[,i]))))

        z.rs =  sapply(1:length(func.labels), function(i) unlist(lm( rank(-log10(pvals.all[,i]))~rank.order)[[1]])[2] )
        z.r  =  sapply(1:length(func.labels), function(i) unlist(lm( rank(-log10(pvals.all[,i]))~rank.effect)[[1]])[2] )
       # z.m  =  sapply(1:length(func.labels), function(i) unlist(lm( (-log10(pvals.all[,i]))~effect.sizes)[[1]])[2] )

	if( PERMUTE_FLAG == TRUE ) {
		cor.permutes.s = list()
		for(k in 1:rr){
			cor.permutes.s[[k]] = sapply(1:length(func.labels), function(i) cor( effect.sizes, -log10(pvals.all[sample(rank.order),i] ), method="s") )
		}
		print("permutes.s done~")
		cor.permutes.p = list()
		for(k in 1:rr){
			cor.permutes.p[[k]] = sapply(1:length(func.labels), function(i) cor( effect.sizes, -log10(pvals.all[sample(rank.order),i] ), method="p") )
		}
	        print("permutes.p done~")
		cor.permutes.r = list()
		for(k in 1:rr){
			cor.permutes.r[[k]] = sapply(1:length(func.labels), function(i) cor( rank.effect, rank(-log10(pvals.all[sample(rank.order),i] ))) )
		}
		print("permutes.r done~")
		cor.permutes.rs = list()
		for(k in 1:rr){
			cor.permutes.rs[[k]] = sapply(1:length(func.labels), function(i) cor( rank.order, rank(-log10(pvals.all[sample(rank.order),i] ))) )
		}
                print("permutes.rs done~")
		z.permutes.rs = list()
		for(k in 1:rr){
			z.permutes.rs[[k]] = sapply(1:length(func.labels), function(i) unlist(lm( rank(-log10(pvals.all[sample(rank.order),i]))~rank.order)[[1]])[2] )
	        }
	        print("permutes.zrs done~")
		z.permutes.r = list()
		for(k in 1:rr){
			z.permutes.r[[k]] = sapply(1:length(func.labels), function(i) unlist(lm( rank(-log10(pvals.all[sample(rank.order),i]))~rank.effect)[[1]])[2] )
		}
	         print("permutes.zr done~")

	#	z.permutes.m = list()
	#	for(k in 1:rr){
	#		z.permutes.m[[k]] = sapply(1:length(func.labels), function(i) unlist(lm( (-log10(pvals.all[sample(rank.order),i]))~effect.sizes)[[1]])[2] )
	
        #	}
	#	 print("permutes.zm done~")

		cor.permutes.list = list(cor.permutes.p, cor.permutes.s, cor.permutes.r, cor.permutes.rs,z.permutes.rs, z.permutes.r)
                names(cor.permutes.list) = labels
	}
	else {
           cor.permutes.list = list()
	}
	cor.list = list(cor.p, cor.s, cor.r, cor.rs,z.rs, z.r)
        names(cor.list) = labels

	return(list( cor.list, cor.permutes.list) )
}


calc_slopes_sub <- function(func.labels,effect.sizes, pvals.all, PERMUTE_FLAG=TRUE, rr=100){
	labels = c("Ranked(spearman corr)")

	rank.order = 1:length(effect.sizes)
	rank.effect = rank(effect.sizes)
	cor.s = sapply(1:length(func.labels), function(i) cor( effect.sizes, -log10(pvals.all[,i] ), method="s") )

	if( PERMUTE_FLAG == TRUE ) {
		cor.permutes.s = list()
		s = list()
		for(k in 1:rr){
                        s[[k]] = sample(rank.order)
			# if( s[[k]] == rank.order) {  }

			cor.permutes.s[[k]] = sapply(1:length(func.labels), function(i) cor( effect.sizes, -log10(pvals.all[s[[k]],i] ), method="s") )
		}
		print("permutes.s done~")

		cor.permutes.list = list(cor.permutes.s)
                names(cor.permutes.list) = labels
	}
	else {
           cor.permutes.list = list()
	}
	cor.list = list(cor.s)
        names(cor.list) = labels

	return(list( cor.list, cor.permutes.list) )
}

calc_FDR <- function(cors,cors.permutes, bbd, bb, nn, pp, rr=100){
	h     = hist(cors, breaks=bb, plot=F)
	h.m   = mean(cors, na.rm=T)
	h.all = cors

	bg     = hist( unlist(cors.permutes), plot=F, breaks=bb)
	bg.all = as.numeric( sapply(1:rr, function(k) cors.permutes[[k]]))

        NN = length(bg.all)
	observed = cbind( (rev(cumsum(rev(h$counts)))), h$mids)
	expected = cbind( nn*(rev(cumsum(rev(bg$counts)))/NN), bg$mids)
	FDR = expected[,1]/observed[,1]

	Pt = expected[ min(which(FDR < pp ) ), 2]
	sig = sum(h.all > Pt, na.rm=T)

#	return( list(Pt, sig))
        return( list(Pt, sig, FDR, h$mids))

}


calc_fdrs_slopes <- function(func.labels,effect.sizes, pvals.all, res, bbd=20, bb=(-1e5:1e5)/1e5 ,pp=0.01){
        cor.list = res[[1]]
        cor.permutes.list = res[[2]]
        nn = length(func.labels)
	rank.order = 1:length(effect.sizes)
	rank.effect = rank(effect.sizes)
        ll = length(cor.list)
        labels = names(cor.list)
        h.c   = lapply(1:ll, function(j) get_density( hist(unlist(cor.list[[j]]) ) ) )
        h     = lapply(1:ll, function(j) get_density( hist(unlist(cor.permutes.list[[j]]) ) ) )
        fdr.c = lapply(1:ll, function(j) calc_FDR( cor.list[[j]], cor.permutes.list[[j]], bbd, bb, nn , pp) )

        par(mfrow=c(2,4))
        for(j in 1:ll){
                ymax = max(h.c[[j]][,2], h[[j]][,2])
                plot(h[[j]], type="l", main=labels[j], xlim=c(-1,1), ylim=c(0,ymax))
                polygon(h[[j]], col=makeTransparent(1))
                polygon(h.c[[j]], col=makeTransparent(2))
                abline(v = fdr.c[[j]][[1]], col=2)
        }
        barplot( unlist(sapply(1:ll, function(j) fdr.c[[j]][[2]])), names=labels , ylab="Significant slopes")
        return(fdr.c)
}



plot_slopes <- function(func.labels,effect.sizes, pvals.all, res){
        nn = length(func.labels)
	rank.order = 1:length(effect.sizes)/length(effect.sizes)
	rank.effect = rank(effect.sizes)/length(effect.sizes)
         pvals.all= as.matrix(pvals.all)
        labels.filt = rownames(pvals.all)

	plot(  0, 0, pch=19, col=0, xlim=c(0,1), ylim=c(0,1), xlab="Ranked disease effect size", ylab="Ranked effect size (sampled n)", axes=F)

	for( i in 1:length(func.labels)) {
		ps = rank(-log10(pvals.all[,i]))/length(rank.effect)
		z = lm( ps~rank.effect)
		points( rank.effect,ps , pch=19 )
		abline(z, lwd=3, col=cols[i])
	}
		axis(1, lab=F)
		text( x=rank.effect, y=0, xpd=T, lab=labels.filt, srt=45)
		axis(2)
}
plot_slopes_raw <- function(func.labels,effect.sizes, pvals.all, res){
        nn = length(func.labels)
	rank.order = 1:length(effect.sizes)/length(effect.sizes)
	rank.effect = rank(effect.sizes)/length(effect.sizes)
	pvals.all= as.matrix(pvals.all)
        labels.filt = rownames(pvals.all)

	plot(  0, 0, col=0, pch=19, xlim=range(effect.sizes), ylim=range(-log10(pvals.all)), xlab="Disease effect size", ylab="Effect size (sampled n)", axes=F)

	for( i in 1:length(func.labels)) {
		ps = -log10(pvals.all[,i])
		z = lm( ps~effect.sizes)
		points( effect.sizes,ps , pch=19 )
		abline(z, lwd=3, col=cols[i])
	}
		axis(1, lab=F)
		text( x=effect.sizes, y=0, xpd=T, lab=labels.filt, srt=45)
		axis(2)
}


