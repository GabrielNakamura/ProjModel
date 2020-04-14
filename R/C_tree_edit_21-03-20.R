library(fishtree)
library(here)

comm<- read.table(here::here("data", "comm_Parana-Paraguai.txt"), header= TRUE)
spp<- colnames(comm)[-c(1,2)]
tree<- fishtree_phylogeny(species = spp, type = "chronogram") #tree downloaded from fishtree of life Chang et al 2019
tree_insert<- tree

###checking the tree#####
spp_data<- 1:length(spp)
names(spp_data)<- spp
insert_spp<- treedata_modif(phy = tree_insert, spp_data, warnings = F)$nc$data_not_tree
species_to_genre<- tree_insert$tip.label[match(sub("_.*", "", insert_spp), 
                                             sub("_.*", "", tree_insert$tip.label)
                                             )[!is.na(match(sub("_.*", "", insert_spp),
                                                            sub("_.*", "", tree_insert$tip.label)
                                                            )
                                                      )
                                               ]
                                         ]  #genre that must be added

species_not_genre<- insert_spp[-unlist(
  unique(
    lapply(species_to_genre, 
           function(x){
             which(sub("_.*", "", insert_spp) == unique(sub("_.*", "", x))
                   )
             }
           )
    )
  )
  ]

##download all species###
families<- c("Loricariidae", "Cichlidae", "Characidae", "Heptapteridae", "Synbranchiformes") 
symbranc_tree<- fishtree_phylogeny(rank= "Synbranchiformes")
spp_allFamilies<- fishtree_phylogeny(species= unlist(lapply(lapply(families, function(x) fishtree_phylogeny(rank= x)), function(i) i$tip.label))) #tree with all species from families in comm
spp_nogenre_allFamilies<- unique(spp_allFamilies$tip.label[match( sub("_.*", "", species_not_genre), 
                                                           sub("_.*", "", spp_allFamilies$tip.label)
                                                           )
                                                    [!is.na(match( sub("_.*", "", species_not_genre), sub("_.*", "", spp_allFamilies$tip.label)))
                                                      ]
                                                    ]
                                 ) #sister taxa on all families
spp_synbranc<- sample(symbranc_tree$tip.label, 1) #extracting one species from Symbranchiformes
include_spp<- c(spp_nogenre_allFamilies, spp_synbranc) #names to be included in spp data
spp_all<- c(spp, include_spp) #all genus to download from FishTree
phylo_raw<- fishtree_phylogeny(species = spp_all) #all genus included 


species_not_genre #position
species_to_genre #add species to genus
spp_problem1<- "Synbranchus_marmoratus"
spp_problem2<- "Pyxiloricaria_menezesi"
pos_1<- treedata_modif(phy = phylo_raw, data = spp_data, warnings = F)$nc$tree_not_data #add position (change name) and species to genus (already had a name in the phylo)
pos_1<- pos_1[-which(pos_1 == spp_synbranc)] #excluding species that are not present in the FishTree 
phylo_test<- force.ultrametric(phylo_raw)

for(i in 1:length(pos_1)){
  #i= 6
  if((table(sub("_.*", "", names(spp_data)))[sub("_.*", "", pos_1[i])] == 1) == TRUE){
    position<- which(sub("_.*", "", phylo_test$tip.label) == sub("_.*", "", pos_1[i]))
    phylo_test$tip.label[position]<- names(spp_data[which(sub("_.*", "", names(spp_data)) == 
                                                           names(table(sub("_.*", "", names(spp_data)))[sub("_.*", "", pos_1[i])] == 1))])
  } else {
    position<- which(sub("_.*", "", phylo_test$tip.label) == sub("_.*", "", pos_1[i]))
    names_spp<- names(spp_data[which(sub("_.*", "", names(spp_data)) == 
                                                           names(table(sub("_.*", "", names(spp_data)))[sub("_.*", "", pos_1[i])] == 1))])
    phylo_test$tip.label[position]<- names_spp[1]
    for(j in 2:length(names_spp)){
      phylo_test<- add.species.to.genus(phylo_test, names_spp[j])
    }
  }  
}

pos_2<- treedata_modif(phy = phylo_test, data = spp_data, warnings = F)$nc$data_not_tree #only species to genus
pos_2<- pos_2[-match(c(spp_problem1, spp_problem2), pos_2)] #symbranchus
for(i in 1:length(pos_2)){
  phylo_test<- add.species.to.genus(phylo_test, pos_2[i])
}


position_problem1<- which(phylo_test$tip.label == spp_synbranc)
phylo_test$tip.label[position_problem1]<- spp_problem1 #solving problem 1

plot(phylo_test, cex= 0.5)
nodelabels()
phylo_test<- bind.tip(tree = phylo_test, tip.label = spp_problem2, where = 75, position = 5)
write.tree(phylo_test, here::here("data", "tree_update_20-03-20.new"))
