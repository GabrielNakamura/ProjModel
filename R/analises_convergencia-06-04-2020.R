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

pcoa.total.taxonomic<- cmdscale(sqrt(total.taxonomic), eig = T) #PCoA in taxonomic turnover matrix
scores.total.taxonomic<- scores(pcoa.total.taxonomic) #extracting scores


#######Phylogenetic beta diversity#######
beta.phylogenetic<-phylo.beta.pair(comm.presau,tree = phylo,index.family="sorensen") #phylogenetic beta diversity
total.phylogenetic<- beta.phylogenetic$phylo.beta.sor
turn.phylogenetic<-beta.phylogenetic$phylo.beta.sim #turnover component of phylogenetic beta diversity
nest.phylogenetic<-beta.phylogenetic$phylo.beta.sne #nestedness component of phylogenetic beta diversity

#phylogenetic total beta div
pcoa.total.phylogenetic<- cmdscale(sqrt(total.phylogenetic), eig = T) #PCoA in phylogenetic turnover matrix
scores.total.phylogenetic<- scores(pcoa.total.phylogenetic)


#phylogenetic turnover
pcoa.turn.phylogenetic<- cmdscale(sqrt(turn.phylogenetic), eig = T) #PCoA in phylogenetic turnover matrix
scores.turn.phylogenetic<- scores(pcoa.turn.phylogenetic)

#phylogenetic nestedness
pcoa.nest.phylogenetic<-cmdscale(sqrt(nest.phylogenetic), eig = T)
scores.nest.phylogenetic<-scores(pcoa.nest.phylogenetic)


#####adonis with environmental factors#####

env_std<- scale(x = env[, -10], center = T, scale = T) #standardizing and centering variables
env_std<- data.frame(cbind(env_std, bacia= as.factor(env$Bacia))) #converting to data frame format


####Adonis testing taxonomic composition####
total_taxonomic_org<- as.matrix(as.dist(total.taxonomic, diag = T, upper = T))[match(rownames(env), rownames(as.matrix(as.dist(total.taxonomic, diag = T, upper = T)))),
                                                                                    match(rownames(env), colnames(as.matrix(as.dist(total.taxonomic, diag = T, upper = T))))]

result_AdonisTaxonomic<- adonis2(total_taxonomic_org~Turbidez + pH + Condutividade + O2.porc + Altitude + Temp + Profundidade + Velocidade + Bacia, 
        data= env, permutations = 999, by = "margin")

result_linear_modTaxonomicBeta<- lm(scores.total.taxonomic[,1] ~ Turbidez + pH + Condutividade + O2.porc + Altitude + Temp + Profundidade + Velocidade + Bacia, 
   data= env)

saveRDS(object = result_AdonisTaxonomic, file = here::here("output", "result_AdonisTaxonomic.rds"))
saveRDS(object = result_linear_modTaxonomicBeta, file = here::here("output", "result_linear_Model_Taxonomic.rds"))


#####Adonis and linear model with taxonomic turnover#####
turn_taxonomic_org<- as.matrix(as.dist(turn.taxonomic, diag = T, upper = T))[match(rownames(env), rownames(as.matrix(as.dist(turn.taxonomic, diag = T, upper = T)))),
                                                                             match(rownames(env), colnames(as.matrix(as.dist(turn.taxonomic, diag = T, upper = T))))]

result_AdonisTaxonomic_turnover<- adonis2(turn_taxonomic_org~Turbidez + pH + Condutividade + O2.porc + Altitude + Temp + Profundidade + Velocidade + Bacia, 
        data= env, permutations = 999, by = "margin") #adonis

result_linear_modTaxonomic_turnover<- lm(scores.total.taxonomic[,1] ~ Turbidez + pH + Condutividade + O2.porc + Altitude + Temp + Profundidade + Velocidade + Bacia, 
                                    data= env)


saveRDS(object = c(result_AdonisTaxonomic, result_linear_modTaxonomic_turnover), file = here::here("output", "result_model_BetaTaxonomic.rds")) #saving specific objects



