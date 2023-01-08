#!/n/app/R/3.5.1/bin/R

args=commandArgs(trailingOnly=TRUE)

library(preprocessCore)
library(RColorBrewer)
library(pheatmap)
library(ggfortify)
library(ggrepel)
library(corrgram)

#percent turned to decimal and negated
perc=1-as.numeric(args[4])/100
#parse hierarchical argument
hier=args[5]
#number of clusters for kmeans, 1 means not
kclusters=as.numeric(args[6])

#A function which takes in a row and returns a 'variability' value
#variability is the difference between EACH value and the mean squared, then averaged
RowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

#function which takes in a row and sets the minimum to 0 and maximum to 1, scales between
RowScale = function(x) {
  (x-min(x))/(max(x)-min(x))
}


#reads in the table properly
a=read.table(args[2],header=TRUE,row.names=1)

#takes the log values of the input table, and runs quantile if necessary
if (as.numeric(args[3])) {
	logvals=log2(as.data.frame(normalize.quantiles(as.matrix(a)))+1)
} else {
	logvals=log2(as.matrix(a)+1)
}


rownames(logvals)=rownames(a)
colnames(logvals)=colnames(a)

#directory for all of the pairwise plots
dir.create(paste(args[1],"_comparison-plots",sep=""), showWarnings = FALSE)

#pairwise dot plot loop, for every single pair x and y... do...
for (x in seq(1,ncol(logvals))){
	for (y in seq(1,ncol(logvals))){
		if (x < y) {
			png(paste(args[1],"_comparison-plots/correlation_",colnames(logvals)[x],"_",colnames(logvals)[y],".png",sep=""))
			#maximum x and y value is set to include at least 99 percent of points, usually more
			#the 1 is an artificial increase
			theMax=max(quantile(logvals[,c(x,y)][,1],0.99),quantile(logvals[,c(x,y)][,2],0.99))+1
			#draws the dot plot using theMax and showing a correlation from the ORIGINAL table, NOT LOG, NOT QUANTILE
			plot(logvals[,c(x,y)],xlim=c(0,theMax),ylim=c(0,theMax),main=cor(a[,c(x,y)])[1,2],cex=0.3)
			dev.off()
		}
	}
}

#the main pdf that all of the plots will go into:
pdf(paste(args[1],"_rplots.pdf",sep=""))

#corrgram shows correlation of samples & upper/lower confidence.  Could be upgraded to adjust size as permitting
corrgram(logvals,lower.panel=panel.conf,upper.panel=panel.cor,cex.labels=0.35,cex=1)

#calculate principal components and HOW MUCH DEVIATION THEY EXPLAIN in PC1 and PC2
prvals=prcomp(as.matrix(t(logvals)))
oneperc=round(prvals$sdev^2/sum(prvals$sdev^2)*100,digits=2)[1]
twoperc=round(prvals$sdev^2/sum(prvals$sdev^2)*100,digits=2)[2]
prvals=prvals$x

#principal component plotting, excess scribbling sets up the boxes and arrows
ggplot(prvals) +
	ggtitle("PCA_all-regions") +
	geom_point(aes(PC1,PC2), size=5,color='grey') +
	geom_label_repel(aes(PC1,PC2,label=colnames(logvals)), fontface='bold',box.padding=unit(0.35,"lines"),point.padding=unit(0.5,"lines"),segment.color='grey50') +
	theme_classic(base_size=16) +
	xlab(paste("PC1 - ",oneperc,sep="")) +
	ylab(paste("PC2 - ",twoperc,sep=""))

#plot a distance matrix, hierarchically cluster or not
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(as.matrix(dist(t(logvals))),col=colors,main="Dist_all",cluster_rows=as.numeric(hier),cluster_cols=as.numeric(hier))

#plot a correlation matrix, hierarchically cluster or not
colors <- colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)
pheatmap(as.matrix(cor(logvals)),col=colors,main="Correlation_all",cluster_rows=as.numeric(hier),cluster_cols=as.numeric(hier))

#number the rows so that original order may be restored later, and then order the table by variability
temp=cbind(logvals,seq(1,nrow(logvals)))[order(RowVar(logvals)),]
#cut the table down to the percent requested
temp=temp[seq(round(nrow(temp)*perc),nrow(temp)),]
#restore the original order of the table, and remove the numbered column
logvals=temp[order(temp[,ncol(temp)]),-ncol(temp)]

#if you DID cut by percent, create set of plots 
#Because PCA already uses the most differential regions... the PCA is probably unnessary and at worst dumb
if (!identical(perc,0)){
	prvals=prcomp(as.matrix(t(logvals)))
	oneperc=round(prvals$sdev^2/sum(prvals$sdev^2)*100,digits=2)[1]
	twoperc=round(prvals$sdev^2/sum(prvals$sdev^2)*100,digits=2)[2]
	prvals=prvals$x


	ggplot(prvals) +
		ggtitle("PCA_variable-regions") +
		geom_point(aes(PC1,PC2), size=5,color='grey') +
		geom_label_repel(
		aes(PC1,PC2,label=colnames(logvals)), fontface='bold',box.padding=unit(0.35,"lines"),point.padding=unit(0.5,"lines"),segment.color='grey50') +
		theme_classic(base_size=16) +
		xlab(paste("PC1 - ",oneperc,sep="")) +
		ylab(paste("PC2 - ",twoperc,sep=""))

	colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
	pheatmap(as.matrix(dist(t(logvals))),col=colors,main="Dist_top-variable",cluster_rows=as.numeric(hier),cluster_cols=as.numeric(hier))

	colors <- colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)
	pheatmap(as.matrix(cor(logvals)),col=colors,main="Correlation_top-variable",cluster_rows=as.numeric(hier),cluster_cols=as.numeric(hier))
}

#if there is to be no clustering, give me some heatmaps
if (identical(kclusters,1)){
		pheatmap(logvals,cluster_rows=F,col=brewer.pal(9,"Blues"),labels_row="",main="Heatmap_top-variable",cluster_cols=as.numeric(hier))
		#row scaled heatmap
		logvalsorder=t(apply(logvals,1,RowScale))
		pheatmap(logvals,cluster_rows=F,col=brewer.pal(9,"Blues"),labels_row="",main="Heatmap_top-variable_row-scaled",cluster_cols=as.numeric(hier))
} else {
	#If you must cluster, here it goes
	#add the cluster number as a column to the table
	temp=cbind(logvals,kmeans(logvals,kclusters)$cluster)
	#extract the order of those clusters
	ordervals=order(temp[,ncol(temp)])
	#actually order by the clusters
	logvalsorder=logvals[ordervals,]
	#set the maximum value of anything to 10 (PROBABLY SHOULD CHANGE OR FIX)
	logvalsorder[logvalsorder>10]=10
	#grab the cluster column in the same ordering function for plotting later
	clusteramounts=as.data.frame(temp[ordervals,ncol(temp)])
	#clusteramounts will be used as annotation for clusters, these step helps color/name them
	colnames(clusteramounts)="cluster"
	rownames(clusteramounts)=rownames(logvalsorder)
	mat_colors=list(cluster = brewer.pal(kclusters,"Paired"))
	names(mat_colors$cluster) = seq(1,kclusters)

	# draw the heatmap with colored/labeled clusters and then again with rowscaling
	pheatmap(logvalsorder,cluster_rows=F,col=brewer.pal(9,"Blues"),annotation_row=clusteramounts,annotation_legend=FALSE,annotation_colors=mat_colors,labels_row="",main="Clustered-Heatmap_top-variable",cluster_cols=as.numeric(hier))
	logvalsorder=t(apply(logvalsorder,1,RowScale))
	pheatmap(logvalsorder,cluster_rows=F,col=brewer.pal(9,"Blues"),annotation_row=clusteramounts,annotation_legend=FALSE,annotation_colors=mat_colors,labels_row="",main="Clustered-Heatmap_top-variable_row-scaled",cluster_cols=as.numeric(hier))


}

dev.off()
