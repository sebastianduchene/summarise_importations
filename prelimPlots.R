library(RColorBrewer)
library(lubridate)

pal <- brewer.pal(n = 3, name = 'Set1')


bangladeshSummary <- read.csv('bangladeshImportationSummary.csv', head = T, sep = ',',
                              stringsAsFactors = F)
bangladeshSummary$firstDate <- decimal_date(ymd(bangladeshSummary$firstDate))
bangladeshSummary$lastDate <- decimal_date(ymd(bangladeshSummary$lastDate))
bangladeshSummary$size <- log10(bangladeshSummary$size)
bangladeshSummary$introDate <- bangladeshSummary$firstDate - bangladeshSummary$detectionLag

dim(bangladeshSummary)
head( bangladeshSummary[, 1:5] )
pdf('introductionsSummaryBangladesh.pdf', useDingbats = F, width = 10, height = 10)
plot(range(2020, bangladeshSummary$lastDate),
     range(bangladeshSummary$size), type = 'n', yaxt = 'n', xaxt = 'n',
     ylab = 'Num. genomes per cluster', xlab = 'Date')
axis(1, at = decimal_date(ymd(c('2020-01-01', '2020-04-01', '2020-08-01', '2020-11-01',
                                '2021-01-01', '2021-04-01', '2021-06-01'))),
     labels = c('Jan \'20', 'Apr \'20',
                'Aug \'20', 'Nov \'20', 'Jan \'21', 'Apr \'21', 'Jun \'21'))
axis(2, at = log10(c(1, 5, 10, 20, 50, 100, 150, 200, 250)),
     labels = c(1, 5, 10, 20, 50, 100, 150, 200, 250))
legend(x = 2020, log10(200), legend = c('first genome', 'last genome', 'infered introduction'),
       bty = 'n', col = c(pal[2:3], 'black'), pch = 20, cex = 1.5)
for(i in 1:nrow(bangladeshSummary)){
    jit <- rnorm(1, mean = 0, sd = 0.05)
    if(bangladeshSummary$size[i] > 0){
        lines(c(bangladeshSummary$introDate[i], bangladeshSummary$lastDate[i]),
              rep(bangladeshSummary$size[i], 2) + jit)
        points(bangladeshSummary$firstDate[i], bangladeshSummary$size[i] + jit,
               col = pal[2], pch = 20, cex = 2)
        points(bangladeshSummary$lastDate[i], bangladeshSummary$size[i] + jit,
               col = pal[3], pch = 20, cex = 2)
        points(bangladeshSummary$introDate[i], bangladeshSummary$size[i] + jit,
               col = 'black', pch = 20, cex = 1)
    }else{
        points(bangladeshSummary$firstDate[i], bangladeshSummary$size[i] + jit,
               col = pal[2], pch = 20, cex = 2)
    }
}
dev.off()

# Summary stats:
source('get_clades.R')
nrow(bangladeshSummary)
sum(bangladeshSummary$size > 0)
date_decimal(range(bangladeshSummary$introDate))
365*mean(bangladeshSummary$detectionLag[bangladeshSummary$size > 0])
get_hindex(10^bangladeshSummary$size)
