library(NELSI)
library(lubridate)
source('get_clades.R')

tr <- read.nexus('nextStrainBangladeshAligned.fasta.treefile.result.date.nexus')
lastTip <- max(as.numeric(gsub('.+_', '', tr$tip.label)))

clades <- find.monophyletic(tr, tag = '_Bangladesh_', include.singletons = T)
nodeDates <- lastTip - intnode.times(tr)
range(nodeDates)

cladesSize <- sapply(clades, function(x) length(x))
sort(cladesSize, decreasing = T)

summaryMatrix <- matrix(NA, length(clades), 6)
colnames(summaryMatrix) <- c('size', 'firstDate', 'lastDate', 'detectionLag', 'lineageCount', 'ids')

for(i in 1:length(clades)){
    dates <- as.numeric(gsub('.+_', '', clades[[i]]))
    detLag <- max(dates) - nodeDates[names(nodeDates) == get.mrca(tr, clades[[i]])]
                                        # Get lineage count and ids
    lineages <- table(gsub('_.+', '', clades[[i]]))
    lineages <- paste(names(lineages), lineages, sep = '=', collapse = ';')
    ids <- paste(clades[[i]], collapse = ';')
    summaryMatrix[i, ] <- c(length(clades[[i]]),
        gsub(' .+', '', c(date_decimal(min(dates)), date_decimal(max(dates)))),
        detLag, lineages, ids)
}
head(summaryMatrix)

write.table(summaryMatrix, file = 'bangladeshImportationSummary.csv', sep = ',', row.names = F)

# Export the top 5 largest clades
aln <- read.dna('nextStrainBangladeshAligned.fasta', format = 'fasta')

top5 <- order(cladesSize, decreasing = T)
for(i in 1:5){
    write.dna(aln[rownames(aln) %in% clades[[ top5[i] ]], ], file = paste0('top_', i, '_clade.fasta'),
              format = 'fasta', nbcol = -1, colsep = '')
}
