library(NELSI)
source("get_clades.R")

#Highest-clade-creibility tree
hcc_tree <- read.nexus('NZ_genomes_V3_SC_CC.tre')

#Some summary data from the tree
youngest_date <- max(as.numeric(gsub('.+[|]', '', hcc_tree$tip.label)))
tips <- hcc_tree$tip.label
node_times <- youngest_date - allnode.times(hcc_tree, reverse = F)

#Transmission lineages and importations are those with the following regular expression tag:
tag <- 'Locally|NA|Import_related|New_Zealand|invest'

# Find monophyletic groups and singletons with the tag above
clades <- find.monophyletic(hcc_tree, tag, include.singletons = T)

# Get closest sister taxon, if the sister is a clade take the oldest tip
first_sisters <- vector()
first_sister_dates <- vector()
for(i in 1:length(clades)){
    sister_temp <- hcc_tree$tip.label[find.sister(hcc_tree, which(hcc_tree$tip.label %in% clades[[i]]))]
    sister_dates <- as.numeric(gsub('.+[|]', '', sister_temp))
    first_sisters[i] <- sister_temp[which.min(sister_dates)]
    first_sister_dates[i] <- min(sister_dates)
}

# Get tmrcas of clades
clades_indexed <- lapply(clades, function(x) which(hcc_tree$tip.label %in% x))
tmrcas <- sapply(clades_indexed,
                 function(x)
                     node_times[get.mrca(hcc_tree, x)])

# Find oldest and youngest sampling times for each clade
first_samples <- sapply(clades_indexed, function(x) min(node_times[x]))
last_samples <- sapply(clades_indexed, function(x) max(node_times[x]))

# Number of samples per clade
num_samples <- sapply(clades, function(x) length(x))

#Make summary matrix
tree_summary <- matrix(NA, length(clades), 7)
colnames(tree_summary) <- c('tmrca', 'first_sample_date',
                            'last_sample_date', 'num_samples',
                            'samples', 'first_sister_tip',
                            'first_sister_date')
tree_summary[, 1] <- tmrcas
tree_summary[, 2] <- first_samples
tree_summary[, 3] <- last_samples
tree_summary[, 4] <- num_samples
tree_summary[, 5] <- sapply(clades, function(x) paste0(x, collapse = ';'))
tree_summary[, 6] <- first_sisters
tree_summary[, 7] <- first_sister_dates
dim(tree_summary)
write.table(tree_summary, file = 'hcc_tree_clades_summary_V3.csv',
            sep = ',', quote = F, row.names = F)



## Now process for 1000 trees from the posterior:
trees <- tail(read.nexus('NZ_genomes_V3_SC_CC.trees'), 1000)

credible_summary <- matrix(NA, length(trees), 3)
colnames(credible_summary) <- c('num_singletons', 'num_lineages', 'mean_size')
for(tr in 1:length(trees)){
    hcc_tree <- trees[[tr]]
                                        #
    youngest_date <- max(as.numeric(gsub('.+[|]', '', hcc_tree$tip.label)))
                                        #
    tips <- hcc_tree$tip.label
    node_times <- youngest_date - allnode.times(hcc_tree, reverse = F)
                                        #
    tag <- 'Locally|NA|Import_related|New_Zealand|invest'
                                        #tag <- 'Locally|NA|Importrelated|New_Zealand|invest|Import'
                                        #
    clades <- find.monophyletic(hcc_tree, tag, include.singletons = T)
                                        #
    first_sisters <- vector()
    first_sister_dates <- vector()
    for(i in 1:length(clades)){
        sister_temp <- hcc_tree$tip.label[find.sister(hcc_tree, which(hcc_tree$tip.label %in% clades[[i]]))]
        sister_dates <- as.numeric(gsub('.+[|]', '', sister_temp))
        first_sisters[i] <- sister_temp[which.min(sister_dates)]
        first_sister_dates[i] <- min(sister_dates)
    }
                                        #
    clades_indexed <- lapply(clades, function(x) which(hcc_tree$tip.label %in% x))
    tmrcas <- sapply(clades_indexed,
                     function(x)
                         node_times[get.mrca(hcc_tree, x)])
    first_samples <- sapply(clades_indexed, function(x) min(node_times[x]))
                                        #
    last_samples <- sapply(clades_indexed, function(x) max(node_times[x]))
    num_samples <- sapply(clades, function(x) length(x))
    tree_summary <- matrix(NA, length(clades), 7)
    colnames(tree_summary) <- c('tmrca', 'first_sample_date',
                                'last_sample_date', 'num_samples',
                                'samples', 'first_sister_tip',
                                'first_sister_date')
    tree_summary[, 1] <- tmrcas
    tree_summary[, 2] <- first_samples
    tree_summary[, 3] <- last_samples
    tree_summary[, 4] <- num_samples
    tree_summary[, 5] <- sapply(clades, function(x) paste0(x, collapse = ';'))
    tree_summary[, 6] <- first_sisters
    tree_summary[, 7] <- first_sister_dates
    credible_summary[tr, ] <- c(sum(tree_summary[, 'num_samples'] == '1'),
                                sum(tree_summary[, 'num_samples'] != '1'),
                                mean(as.numeric(tree_summary[, 'num_samples'])))
    print(tr)
    print(credible_summary[tr, ])
}

write.table(credible_summary, file = 'credible_summary_V3.csv', row.names = F,
            sep = ',', quote = F)
