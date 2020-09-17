library(NELSI)
library(lubridate)
library(RColorBrewer)

post_file <- tail(read.table('BD_paper.afa.datedR2_estimate_cut.log',  head = T), 1500)
tr <- read.nexus('BD_paper.afa.datedR2_estimate_cut.tre')

r1 <- density(post_file$reproductiveNumber_BDSKY_Serial.1, bw = 5e-2)
r2 <- density(post_file$reproductiveNumber_BDSKY_Serial.2, bw = 5e-2)
ct <- density(post_file$ReChangeTime2, bw = 1e-3)


th <- max(allnode.times(tr))
num_tips <- length(tr$tip.label)
latest_tip <- as.numeric(max(gsub('.+@', '', tr$tip.label)))


median_cut_height <- th -  quantile(ct$x, 0.5)
tip_heights <- allnode.times(tr, tipsonly = T)
tips_pre <- which(tip_heights <= median_cut_height)
tips_post <- which(tip_heights > median_cut_height)

blues <- brewer.pal(9, 'Blues')
oranges <- brewer.pal(9, 'Oranges')

#pdf('paper-20200424/BD_sky_estimate_Rechange_R2_v2.pdf', useDingbats = F,
#    height = 5, width = 9)
par(mar = c(2, 2, 2, 2))
plot(ladderize(tr), show.tip.label = F, edge.col = 'darkgrey')
tiplabels(text = rep('', length(tips_pre)),
          tip = tips_pre, frame = 'none', cex = 0.5, pch = 20,
          col = blues[5])
tiplabels(text = rep('', length(tips_pre)),
          tip = tips_post, frame = 'none', cex = 0.5, pch = 20,
          col = blues[9])
blue1 <- col2rgb(blues[7], alpha = T)/255
polygon(th-ct$x, ct$y*4, col = rgb(blue1[1], blue1[2], blue1[3], 0.6),
        border =  blues[7], lwd = 2)
polygon((r1$y/200)+th/2, r1$x*450, col = rgb(1, 1, 1, 0.5),
        border = 'black', lwd = 2)
polygon(-(r1$y/200)+th/2, r1$x*450, col = rgb(1, 1, 1, 0.5),
        border = 'black', lwd = 2)
points(th/2, r1$x[which.max(r1$y)]*450, pch = 3, lwd = 2)
polygon((r2$y/200)+th/1.05, r2$x*450, col = rgb(1, 1, 1, 0.5),
        border = 'black', lwd = 2)
polygon(-(r2$y/200)+th/1.05, r2$x*450, col = rgb(1, 1, 1, 0.5),
        border = 'black', lwd = 2)
points(th/1.05, r2$x[which.max(r2$y)]*450, pch = 3, lwd = 2)

