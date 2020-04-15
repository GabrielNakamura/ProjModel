######ordination graphics

####total taxonomic beta diversity#####
quartz() #open a new window to plot ordination graphic
plot(scores.total.taxonomic[,1],scores.total.taxonomic[,2],type="n",
     main="PCoA taxonomic beta-div")
points(scores.total.taxonomic[1:10,1],scores.total.taxonomic[1:10,2],pch=19,cex=1.8)
points(scores.total.taxonomic[11:20,1],scores.total.taxonomic[11:20,2],pch=2,cex=1.8)
text(scores.total.taxonomic[,1],scores.total.taxonomic[,2], labels=rownames(scores.total.taxonomic))
abline(v=0,h=0,lty=2)

pdf(file = here::here("output", "figures", "Fig01_taxonomic_turnover.pdf"))
plot(scores.total.taxonomic[,1],scores.total.taxonomic[,2],type="n",
     main="PCoA taxonomic beta-div")
points(scores.total.taxonomic[1:10,1],scores.total.taxonomic[1:10,2],pch=19,cex=1.8)
points(scores.total.taxonomic[11:20,1],scores.total.taxonomic[11:20,2],pch=2,cex=1.8)
text(scores.total.taxonomic[,1],scores.total.taxonomic[,2], labels=rownames(scores.total.taxonomic))
abline(v=0,h=0,lty=2)
dev.off()

#####Taxonomic turnover#######
pcoa.turn.taxonomic<- cmdscale(sqrt(turn.taxonomic), eig = T) #PCoA in taxonomic turnover matrix
scores.turn.taxonomic<- scores(pcoa.turn.taxonomic)
quartz()
plot(scores.turn.taxonomic[,1],scores.turn.taxonomic[,2],type="n",
     main="PCoA taxonomic turnover")
points(scores.turn.taxonomic[1:10,1],scores.turn.taxonomic[1:10,2],pch=19,cex=1.8)
points(scores.turn.taxonomic[11:20,1],scores.turn.taxonomic[11:20,2],pch=2,cex=1.8)
text(scores.turn.taxonomic[,1],scores.turn.taxonomic[,2], labels=rownames(scores.turn.taxonomic))
abline(v=0,h=0,lty=2)



#####Phylogenetic beta diversity######
quartz()
plot(scores.total.phylogenetic[,1],scores.total.phylogenetic[,2],type="n",
     main="PCoA phylogenetic beta-div")
points(scores.total.phylogenetic[1:10,1],scores.total.phylogenetic[1:10,2],pch=19,cex=1.8)
points(scores.total.phylogenetic[11:20,1],scores.total.phylogenetic[11:20,2],pch=2,cex=1.8)
text(scores.total.phylogenetic[,1], (scores.total.phylogenetic[,2] + 0.02), labels=rownames(scores.total.phylogenetic))
abline(v=0,h=0,lty=2)

#####Phylogenetic turnover######
pcoa.turn.phylogenetic<- cmdscale(sqrt(turn.phylogenetic), eig = T) #PCoA in phylogenetic turnover matrix
scores.turn.phylogenetic<- scores(pcoa.turn.phylogenetic)
quartz()
plot(scores.turn.phylogenetic[,1],scores.turn.phylogenetic[,2],type="n",
     main="PCoA phylogenetic turnover")
points(scores.turn.phylogenetic[1:10,1],scores.turn.phylogenetic[1:10,2],pch=19,cex=1.8)
points(scores.turn.phylogenetic[11:20,1],scores.turn.phylogenetic[11:20,2],pch=2,cex=1.8)
text(scores.turn.phylogenetic[,1], (scores.turn.phylogenetic[,2] + 0.02), labels=rownames(scores.turn.phylogenetic))
abline(v=0,h=0,lty=2)