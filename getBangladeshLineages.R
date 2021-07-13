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

summaryMatrix <- matrix(NA, length(clades), 5)
colnames(summaryMatrix) <- c('firstDate', 'lastDate', 'detectionLag', 'lineageCount', 'ids')

for(i in 1:length(clades)){
    dates <- as.numeric(gsub('.+_', '', clades[[i]]))
    detLag <- max(dates) - nodeDates[names(nodeDates) == get.mrca(tr, clades[[i]])]
    # Get lineage count and ids
    summaryMatrix[i, 1:3] <- c(gsub(' .+', '',
                                    c(date_decimal(min(dates)), date_decimal(max(dates)))),
                              detLag)
}
head(summaryMatrix)
