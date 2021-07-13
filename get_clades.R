library(NELSI)
library(geiger)


get_hindex <- function(x){
#    x <- sample(1:10, size = 10, replace = T)
    sorted <- sort(x, decreasing = T)
    i <- 1:length(sorted)
    hs <- sorted[which(sorted <= i)]
    return(hs[1])
}


get.mrca <- function (phy, tip) 
{
    if (!inherits(phy, "phylo")) 
        stop("object \"phy\" is not of class \"phylo\"")
    if (length(tip) < 2){
        tnd <- if (is.character(tip)) 
                   match(tip, phy$tip.label)
               else tip
        return(phy$edge[phy$edge[, 2] == tnd, 1])
    }
    Ntip <- length(phy$tip.label)
    rootnd <- Ntip + 1L
    pars <- integer(phy$Nnode)
    tnd <- if (is.character(tip)) 
        match(tip, phy$tip.label)
    else tip
    done_v <- logical(Ntip + phy$Nnode)
    pvec <- integer(Ntip + phy$Nnode)
    pvec[phy$edge[, 2]] <- phy$edge[, 1]
    nd <- tnd[1]
    for (k in 1:phy$Nnode) {
        nd <- pvec[nd]
        pars[k] <- nd
        if (nd == rootnd) 
            break
    }
    pars <- pars[1:k]
    mrcind <- integer(max(pars))
    mrcind[pars] <- 1:k
    mrcand <- pars[1]
    for (i in 2:length(tnd)) {
        cnd <- tnd[i]
        done <- done_v[cnd]
        while (!done) {
            done_v[cnd] <- TRUE
            cpar <- pvec[cnd]
            done <- done_v[cpar]
            if (cpar %in% pars) {
                if (cpar == rootnd) 
                  return(rootnd)
                if (mrcind[cpar] > mrcind[mrcand]) 
                  mrcand <- cpar
                done_v[cpar] <- TRUE
                done <- TRUE
            }
            cnd <- cpar
        }
    }
    mrcand
}


get.descendants.of.node <-
function(phy, node, tips=FALSE){
	n=Ntip(phy)
	all=ifelse(tips, FALSE, TRUE)
	out <- .Call("get_descendants", tree=list(
    NODE = as.integer(node),
    ROOT = as.integer(n+1),
    ALL = as.integer(all),
    ENDOFCLADE = as.integer(dim(phy$edge)[1]),
    ANC = as.integer(phy$edge[,1]),
    DES = as.integer(phy$edge[,2])),
    PACKAGE = "geiger")
	res=out$TIPS
	if(!length(res)) res=NULL
	return(c(node, res))
}

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
    require(geiger)
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
#        descending_nodes <- get.descending.nodes.branches(tr, intnodes[i])$descending.nodes
        descending_nodes <- get.descendants.of.node(tr, intnodes[i])
        descending_tips <- tr$tip.label[descending_nodes[descending_nodes %in% tips]]
##        descending_tips <- tips(tr, intnodes[i])
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





