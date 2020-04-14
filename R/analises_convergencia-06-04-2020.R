###MS Taxonomic and phylogenetic turnover in Parana and Paraguay streams

####packages####
library(vegan)
library(ape)
library(betapart)
library(SYNCSA)
library(ade4)
library(geiger)
library(phytools)

####functions####
source(here::here("R", "function", "function_superTree_fish.R"))
result_linear_model_taxonomic<- result_linear_model_taxonomic<- readRDS(file = here::here("output", "result_linear_Model_Taxonomic.rds")) #reading results from linear model with taxonomic beta diversity
result_modelTaxonomic<- readRDS(file = here::here("output", "result_model_BetaTaxonomic.rds")) #reading models (Adonis and lm)


#####data####
comm<- read.table(here::here("data", "processed", "comm_Parana-Paraguai.txt"), header= TRUE) #community data
phylo<- read.tree(here::here("data", "processed", "tree_update_20-03-20.new")) #phylogeny
env<- read.table(here::here("data", "processed", "amb.txt"), header= T) #environmental variables
comm.presau<- ifelse(comm[,3:ncol(comm)]>=1,1,0) #incidence matrix
rownames(comm.presau)<- comm[,"pontos"] #naming the rows with the names of streams

###saving proccessed files####
write.table(comm.presau, here::here("data", "processed", "comm_presau.txt")) #write community occurence data in txt format

#######Taxonomic beta diversity#######
#between basins
beta.taxonomic<-beta.pair(comm.presau,index.family="sorensen") #taxonomic beta diversity
total.taxonomic<- beta.taxonomic$beta.sor #total beta taxonomic 
turn.taxonomic<-beta.taxonomic$beta.sim #turnover component of taxonomic beta diversity
nest.taxonomic<-beta.taxonomic$beta.sne #nestedness component of taxonomic beta diversity

#plotting total beta-div
pcoa.total.taxonomic<- cmdscale(sqrt(total.taxonomic), eig = T) #PCoA in taxonomic turnover matrix
scores.total.taxonomic<- scores(pcoa.total.taxonomic) #extracting scores
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

#plotting turnover
pcoa.turn.taxonomic<- cmdscale(sqrt(turn.taxonomic), eig = T) #PCoA in taxonomic turnover matrix
scores.turn.taxonomic<- scores(pcoa.turn.taxonomic)
quartz()
plot(scores.turn.taxonomic[,1],scores.turn.taxonomic[,2],type="n",
     main="PCoA taxonomic turnover")
points(scores.turn.taxonomic[1:10,1],scores.turn.taxonomic[1:10,2],pch=19,cex=1.8)
points(scores.turn.taxonomic[11:20,1],scores.turn.taxonomic[11:20,2],pch=2,cex=1.8)
text(scores.turn.taxonomic[,1],scores.turn.taxonomic[,2], labels=rownames(scores.turn.taxonomic))
abline(v=0,h=0,lty=2)



#plotting nestedness
pcoa.nest.taxonomic<- cmdscale(sqrt(nest.taxonomic), eig = T)
scores.nest.taxonomic<-scores(pcoa.nest.taxonomic)
quartz()
plot(scores.nest.taxonomic[,1],scores.nest.taxonomic[,2],type="n",
     main="PCoA taxonomic nestedness")
points(scores.nest.taxonomic[1:10,1],scores.nest.taxonomic[1:10,2],pch=19,cex=1.8)
points(scores.nest.taxonomic[11:20,1],scores.nest.taxonomic[11:20,2],pch=2,cex=1.8)
text(scores.nest.taxonomic[,1],scores.nest.taxonomic[,2], labels=rownames(scores.nest.taxonomic))
abline(v=0,h=0,lty=2)

#plot with betadisper function
quartz()
plot(betadisper(d = total.taxonomic, group = comm$bacia, sqrt.dist = T), cex= 1.3, main= "Total taxonomic diversity", 
     xlab = "PCoA 1 (21.8 %)", ylab= "PCoA 2 (12 %)")
permutest(betadisper(d = total.taxonomic, group = comm$bacia, sqrt.dist = T))

quartz()
plot(betadisper(d = turn.taxonomic, group = comm$bacia, sqrt.dist = T), cex= 1.3, main= "Taxonomic turnover", 
     xlab = "PCoA 1 (35.4 %)", ylab= "PCoA 2 (15.9 %)")

quartz()
plot(betadisper(d = nest.taxonomic, group = comm$bacia, sqrt.dist = T), cex= 1.3, main= "Taxonomic nestedness", 
     xlab = "PCoA 1 (66 %)", ylab= "PCoA 2 (49 %)")
permutest(betadisper(d = nest.taxonomic, group = comm$bacia, sqrt.dist = T))



#######Phylogenetic beta diversity#######
beta.phylogenetic<-phylo.beta.pair(comm.presau,tree = phylo,index.family="sorensen") #phylogenetic beta diversity
total.phylogenetic<- beta.phylogenetic$phylo.beta.sor
turn.phylogenetic<-beta.phylogenetic$phylo.beta.sim #turnover component of phylogenetic beta diversity
nest.phylogenetic<-beta.phylogenetic$phylo.beta.sne #nestedness component of phylogenetic beta diversity

#phylogenetic total beta div
pcoa.total.phylogenetic<- cmdscale(sqrt(total.phylogenetic), eig = T) #PCoA in phylogenetic turnover matrix
scores.total.phylogenetic<- scores(pcoa.total.phylogenetic)
quartz()
plot(scores.total.phylogenetic[,1],scores.total.phylogenetic[,2],type="n",
     main="PCoA phylogenetic beta-div")
points(scores.total.phylogenetic[1:10,1],scores.total.phylogenetic[1:10,2],pch=19,cex=1.8)
points(scores.total.phylogenetic[11:20,1],scores.total.phylogenetic[11:20,2],pch=2,cex=1.8)
text(scores.total.phylogenetic[,1], (scores.total.phylogenetic[,2] + 0.02), labels=rownames(scores.total.phylogenetic))
abline(v=0,h=0,lty=2)

#phylogenetic turnover
pcoa.turn.phylogenetic<- cmdscale(sqrt(turn.phylogenetic), eig = T) #PCoA in phylogenetic turnover matrix
scores.turn.phylogenetic<- scores(pcoa.turn.phylogenetic)
quartz()
plot(scores.turn.phylogenetic[,1],scores.turn.phylogenetic[,2],type="n",
     main="PCoA phylogenetic turnover")
points(scores.turn.phylogenetic[1:10,1],scores.turn.phylogenetic[1:10,2],pch=19,cex=1.8)
points(scores.turn.phylogenetic[11:20,1],scores.turn.phylogenetic[11:20,2],pch=2,cex=1.8)
text(scores.turn.phylogenetic[,1], (scores.turn.phylogenetic[,2] + 0.02), labels=rownames(scores.turn.phylogenetic))
abline(v=0,h=0,lty=2)

#phylogenetic nestedness
pcoa.nest.phylogenetic<-cmdscale(sqrt(nest.phylogenetic), eig = T)
scores.nest.phylogenetic<-scores(pcoa.nest.phylogenetic)
quartz()
plot(scores.nest.phylogenetic[,1],scores.nest.phylogenetic[,2],type="n",
     main="PCoA phylogenetic nestedness")
points(scores.nest.phylogenetic[1:10,1],scores.nest.phylogenetic[1:10,2],pch=19,cex=1.8)
points(scores.nest.phylogenetic[11:20,1],scores.nest.phylogenetic[11:20,2],pch=2,cex=1.8)
text(scores.nest.phylogenetic[,1],scores.nest.phylogenetic[,2], labels=rownames(scores.nest.phylogenetic))
abline(v=0,h=0,lty=2)


#####adonis with environmental factors#####

env_std<- scale(x = env[, -10], center = T, scale = T)
env_std<- data.frame(cbind(env_std, bacia= as.factor(env$Bacia)))


#testing taxonomic composition
total_taxonomic_org<- as.matrix(as.dist(total.taxonomic, diag = T, upper = T))[match(rownames(env), rownames(as.matrix(as.dist(total.taxonomic, diag = T, upper = T)))),
                                                                                    match(rownames(env), colnames(as.matrix(as.dist(total.taxonomic, diag = T, upper = T))))]

result_AdonisTaxonomic<- adonis2(total_taxonomic_org~Turbidez + pH + Condutividade + O2.porc + Altitude + Temp + Profundidade + Velocidade + Bacia, 
        data= env, permutations = 999, by = "margin")

result_linear_modTaxonomicBeta<- lm(scores.total.taxonomic[,1] ~ Turbidez + pH + Condutividade + O2.porc + Altitude + Temp + Profundidade + Velocidade + Bacia, 
   data= env)

saveRDS(object = result_AdonisTaxonomic, file = here::here("output", "result_AdonisTaxonomic.rds"))
saveRDS(object = result_linear_modTaxonomicBeta, file = here::here("output", "result_linear_Model_Taxonomic.rds"))


#adonis and linear model with turnover
turn_taxonomic_org<- as.matrix(as.dist(turn.taxonomic, diag = T, upper = T))[match(rownames(env), rownames(as.matrix(as.dist(turn.taxonomic, diag = T, upper = T)))),
                                                                             match(rownames(env), colnames(as.matrix(as.dist(turn.taxonomic, diag = T, upper = T))))]

result_AdonisTaxonomic_turnover<- adonis2(turn_taxonomic_org~Turbidez + pH + Condutividade + O2.porc + Altitude + Temp + Profundidade + Velocidade + Bacia, 
        data= env, permutations = 999, by = "margin") #adonis

result_linear_modTaxonomic_turnover<- lm(scores.total.taxonomic[,1] ~ Turbidez + pH + Condutividade + O2.porc + Altitude + Temp + Profundidade + Velocidade + Bacia, 
                                    data= env)


saveRDS(object = c(result_AdonisTaxonomic, result_linear_modTaxonomic_turnover), file = here::here("output", "result_model_BetaTaxonomic.rds")) #saving specific objects

