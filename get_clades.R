library(NELSI)


find.sister <- function(tr, clade){
    if(!is.monophyletic(tr, clade)) stop('group is not monophyletic')
    mrca <- get.mrca(tr, clade)
    tips <- 1:length(tr$tip.label)
    parent_node <- tr$edge[tr$edge[, 2] == mrca, 1]
    all_descendants <- get.descending.nodes.branches(tr, parent_node)$descending.nodes
    clade_nodes <- get.descending.nodes.branches(tr, mrca)$descending.nodes
    sister_nodes <- all_descendants[!all_descendants %in% clade_nodes]
    sister_tips <- sister_nodes[sister_nodes %in% tips]
    return(sister_tips)
}


is.polytomy <- function(tr, node){
    sum(tr$edge[, 1] == node) > 2
}

find.monophyletic <- function(tr, tag, include.singletons = F){
#include.singletons = F
#tr <- read.tree(text = '(((a, b_uno, c_uno), (d, e, (f, g))), h);')
#tr <- read.tree(text = '((((a, i_uno), b_uno, c_uno), (d, e, (f, g))), h);')
#tr <- read.tree(text = '(((((a, (j_uno, k_uno)), i_uno), b_uno, c_uno), (d, e, (f, g))), h);')
#tr <- read.tree(text = '(((((a, (j_uno, k_uno)), i_uno), b_uno, c_uno), (d, e, (f, g))), (h_uno, l_uno));')
#tag <- 'uno'
#plot(tr)
#tiplabels(cex = 0.7)
#nodelabels()
    if(!is.rooted(tr)) stop('tree is not rooted')
    tips <- 1:length(tr$tip.label)
    intnodes <- unique(tr$edge[, 1])
    clades <- list()
    i <- 1
    while(i <= length(intnodes)){
        descending_nodes <- get.descending.nodes.branches(tr, intnodes[i])$descending.nodes
        descending_tips <- tr$tip.label[descending_nodes[descending_nodes %in% tips]]
        num_matches <- grep(tag, descending_tips)
        num_internal_nodes <- length(descending_nodes) - length(descending_tips)
        if(length(num_matches) == 0){
            intnodes <- intnodes[-which(intnodes %in% descending_nodes)]
        }else if(length(num_matches) == length(descending_tips)){
            clades[[length(clades)+1]] <- descending_tips
            intnodes <- intnodes[-which(intnodes %in% descending_nodes)]
        }else if(is.polytomy(tr, intnodes[i]) ){ #If it is a polytomy
            i_matches <- which(tr$tip.label %in% descending_tips[num_matches])
            ancestral_nodes <- sapply(i_matches, function(x) get.mrca(tr, x))
            tips_in_multi_clade <- i_matches[ancestral_nodes == intnodes[i]]
            if(length(tips_in_multi_clade) > 1){
                clades[[length(clades)+1]] <- tr$tip.label[tips_in_multi_clade]
            }
            intnodes <- intnodes[intnodes != intnodes[i]]
        }else{
            i <- i+1
        }
    }




if(include.singletons){
        tips_in_clades <- unlist(clades)
        tips_out_clades <- tr$tip.label[!tr$tip.label %in% tips_in_clades]
        tips_out_clades <- grep(tag, tips_out_clades, value = T)
        return(c(clades, tips_out_clades))
    }else{
        return(clades)
    }
}



#tr <- read.tree(text = '(((a, b_uno, c_uno), (d, e_uno, (f, g))), h);')
#plot(tr)

#find.monophyletic(tr, '_uno')




#example:
#par(mfrow = c(1, 1))
#set.seed(144)
#tr <- rtree(15)
#tr$tip.label[c(1, 2, 3, 6, 7, 15, 13)] <- paste0(tr$tip.label[c(1, 2, 3, 6, 7, 15, 13)], '_local')
#plot(tr)
#nodelabels()
#tag <- 'local'
#find.monophyletic(tr, 'local', T)







#}

###########
get.clades <- function(tr, threshold){
    'threshold is in forwards in time with the age of the root set to zero'
                                        #    set.seed(124)
                                        #    tr <- rtree(10)
                                        #    threshold <- 1.8
    nts <- allnode.times(tr, reverse = T)
    nts <- nts[!names(nts) %in% as.character(1:length(tr$tip.label))]
                                        #    plot(tr)
                                        #    lines(rep(threshold, 2), c(0, 10))
                                        #     nodelabels(round(nts, 2))
    nodes_past_threshold <- names(nts[nts > threshold])
#    nodelabels(text = nodes_past_threshold,
#               node = as.numeric(nodes_past_threshold),
#               bg = 'green')

    descendants_list <- list()
    for(x in nodes_past_threshold){
        descendants_list[[length(descendants_list) +1]] <-
            get.descending.nodes.branches(tr, x)$descending.nodes[-1]
    }
    names(descendants_list) <- nodes_past_threshold
    contained <- vector()
    clist <- list()
    query_clades <- descendants_list
    for(i in 1:length(descendants_list)){
        for(j in 1:length(query_clades)){
            check_contained <- all(descendants_list[[i]] %in% query_clades[[j]])
            if(check_contained){
                contained <- c(contained, names(query_clades)[j])
                break
            }
        }
    }
    contained <- unique(contained)
    list_taxa <- lapply(contained,
                        function(x)
                            as.numeric(get.descending.nodes.branches(tr, x)$descending.nodes))
    names(list_taxa) <- contained
    return(list_taxa)
}
