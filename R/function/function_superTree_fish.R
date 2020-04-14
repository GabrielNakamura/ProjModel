superTree_fish<- function(comm_data, families, spp_problem1= NULL, spp_problem2= NULL, rank_problem1= NULL, rank_problem2= NULL, where, position){
  
  ####internal function#####
  treedata_modif<- function (phy, data, sort = FALSE, warnings = TRUE) 
  {
    dm = length(dim(data))
    if (is.vector(data)) {
      data <- as.matrix(data)
    }
    if (is.factor(data)) {
      data <- as.matrix(data)
    }
    if (is.array(data) & length(dim(data)) == 1) {
      data <- as.matrix(data)
    }
    if (is.null(rownames(data))) {
      stop("names for 'data' must be supplied")
    }
    else {
      data.names <- rownames(data)
    }
    nc <- name.check(phy, data)
    if (is.na(nc[[1]][1]) | nc[[1]][1] != "OK") {
      if (length(nc[[1]] != 0)) {
        phy = drop.tip(phy, as.character(nc[[1]]))
        if (warnings) {
          warning(paste("The following tips were not found in 'data' and were dropped from 'phy':\n\t", 
                        paste(nc[[1]], collapse = "\n\t"), sep = ""))
        }
      }
      if (length(nc[[2]] != 0)) {
        m <- match(data.names, nc[[2]])
        data = as.matrix(data[is.na(m), ])
        data.names <- data.names[is.na(m)]
        if (warnings) {
          warning(paste("The following tips were not found in 'phy' and were dropped from 'data':\n\t", 
                        paste(nc[[2]], collapse = "\n\t"), sep = ""))
        }
      }
    }
    order <- match(data.names, phy$tip.label)
    rownames(data) <- phy$tip.label[order]
    if (sort) {
      index <- match(phy$tip.label, rownames(data))
      data <- as.matrix(data[index, ])
    }
    if (dm == 2) {
      data <- as.matrix(data)
    }
    phy$node.label = NULL
    return(list(phy = phy, data = data, nc= nc))
  }
  #########
  
  spp<- colnames(comm_data)
  tree_insert<- fishtree_phylogeny(species = spp, type = "chronogram")
  
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
  ] #species without any genre in the tree
  
  ##downloading all species from FishTree
  spp_allFamilies<- fishtree_phylogeny(species= unlist(lapply(lapply(families, 
                                                                     function(x) fishtree_phylogeny(rank= x)), 
                                                              function(i) i$tip.label)
                                                       )
                                       ) #tree with all species from families in comm

  if(!is.null(spp_problem1)){
    families<- c(families, rank_problem1)
    spp_allFamilies<- fishtree_phylogeny(species= unlist(lapply(lapply(families, 
                                                                       function(x) fishtree_phylogeny(rank= x)), 
                                                                function(i) i$tip.label)
    )
    )
  }
  if(!is.null(spp_problem2)){
    families<- c(families, rank_problem2)
    spp_allFamilies<- fishtree_phylogeny(species= unlist(lapply(lapply(families, 
                                                                       function(x) fishtree_phylogeny(rank= x)), 
                                                                function(i) i$tip.label)
    )
    )
  }
  
  spp_nogenre_allFamilies<- unique(spp_allFamilies$tip.label[match( sub("_.*", "", species_not_genre), 
                                                                    sub("_.*", "", spp_allFamilies$tip.label))[
                                                                      !is.na(match( sub("_.*", "", species_not_genre), 
                                                                                    sub("_.*", "", spp_allFamilies$tip.label)
                                                                                    )
                                                                             )
                                                                      ]
                                                             ]
                                   )  #downloading all species from families
  
  if(is.null(spp_problem1)){ #no species on problem 1
    phylo_raw<- fishtree_phylogeny(species = spp_nogenre_allFamilies) #all genus included 
  }
  
  if(!is.null(spp_problem1)){
    if(length(rank_problem1) == 1){
      problem1_tree<- fishtree_phylogeny(rank= rank_problem1)
      spp_samp_problem1<- sample(problem1_tree$tip.label, 1) #extracting one species to solve problem 1
      include_spp<- c(spp_nogenre_allFamilies, spp_samp_problem1) #species name to be included in spp data
      spp_all<- c(spp, include_spp) #all genus to download from FishTree
      phylo_raw<- fishtree_phylogeny(species = spp_all) #all genus included 
    }
  } else{
    spp_samp_problem1<- numeric(length = length(rank_problem1))
    for(i in 1:rank_problem1){
      problem1_tree<- fishtree_phylogeny(rank= rank_problem1[i])
      spp_samp_problem1[i]<- sample(problem1_tree$tip.label, 1)
    }
    include_spp<- c(spp_nogenre_allFamilies, spp_samp_problem1) #species name to be included in spp data
    spp_all<- c(spp, include_spp) #all genus to download from FishTree
    phylo_raw<- fishtree_phylogeny(species = spp_all) #all genus included 
  } #download phylogeny from FishTree to be edited
  
  #selecting species to add and to substitute
  pos_1<- treedata_modif(phy = phylo_raw, data = spp_data, warnings = F)$nc$tree_not_data #add position (change name) and species to genus (already had a name in the phylo)
  if(length(spp_problem1) != 0){
    pos_1<- pos_1[-which(pos_1 == spp_problem1)] #excluding species that are not present in the FishTree 
  } else{
    pos_1 
  }
  
  phylo_test<- force.ultrametric(phylo_raw) #ultrameterizing the phylogeny that will be manipulated
  
  #replacing the names of species from FishTree for species in pool and/or including sister species
  for(i in 1:length(pos_1)){
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
  
  #accessing only species that must be added to genus already existant in phylogenetic tree
  pos_2<- treedata_modif(phy = phylo_test, data = spp_data, warnings = F)$nc$data_not_tree #only species to genus
  
  #removing species in pos2 accordingly problem 1 and problem 2
  if(!is.null(spp_problem1)){
    pos_2<- pos_2[-match(spp_problem1, pos_2)] #removing species from problem 1
  }
  if(!is.null(spp_problem2)){
    pos_2<- pos_2[-match(spp_problem2, pos_2)] #removing species from problem 1
  }
  
  for(i in 1:length(pos_2)){
    phylo_test<- add.species.to.genus(phylo_test, pos_2[i]) #adding species to genus, except the problems 1 and 2
  }
  
  
  if(!is.null(problem1)){
    if(length(problem1) == 1){
      position_problem1<- which(phylo_test$tip.label == spp_samp_problem1)
      phylo_test$tip.label[position_problem1]<- spp_problem1 #solving problem 1
    }
    if(length(spp_problem1 > 1)){ #case in which spp_problem1 have more than one species
      for(i in 2:length(spp_problem1)){
        phylo_test<- add.species.to.genus(phylo_test, spp_samp_problem1)
        position_problem1<- which(phylo_test$tip.label == spp_samp_problem1)
        phylo_test$tip.label[position_problem1]<- spp_problem1[i]
      }
    }
  }
  
  
  #solving problem 2
  if(!is.null(spp_problem2)){
    if(length(spp_problem2) == 1){
      phylo_test<- bind.tip(tree = phylo_test, tip.label = spp_problem2, where = where, position = position) #solving for problem 2
      if(length(spp_problem2 > 1)){ #solving problem 2 for more than one species
        for( i in 1:length(spp_problem2)){
          phylo_test<- bind.tip(tree = phylo_test, tip.label = spp_problem2[i], where = where[i], position = position[i]) #solving for problem 2
        }
      }
    }
  }

  
  phylo_test #supertree
  
}
