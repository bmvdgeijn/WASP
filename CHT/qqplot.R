



filenames <- c(# "/iblm/netapp/home/gmcvicker/proj/WASP/CHT/nick_input_data_chr22/cht_results.txt",
               "/iblm/netapp/home/gmcvicker/proj/WASP/CHT/nick_input_data_chr22/cht_results_as.txt",
               "/iblm/netapp/home/gmcvicker/proj/WASP/CHT/nick_input_data_chr22/cht_results_bnb.txt",
              # "/iblm/netapp/home/gmcvicker/proj/WASP/CHT/nick_input_data_chr22/cht_results_permuted.txt",
                "/iblm/netapp/home/gmcvicker/proj/WASP/CHT/nick_input_data_chr22/cht_results_as_permuted.txt",
               "/iblm/netapp/home/gmcvicker/proj/WASP/CHT/nick_input_data_chr22/cht_results_bnb_permuted.txt"
               )


library(RColorBrewer)

pal <- brewer.pal(9, "Set1")

png("qqplot.png", width=500, height=500)

labels <- c()

min.p <- 1e-20

for(i in 1:length(filenames)) {
    filename <- filenames[i]

    tab <- read.table(filename, header=T)
    n.test <- nrow(tab)
    null.p <- (1:n.test)/(n.test)
    obs.p <- tab$P.VALUE

    # cap min p for drawing purposes
    obs.p[obs.p < min.p]  <- min.p
    null.p[null.p < min.p] <- min.p

    vals <- qqplot(-log10(null.p), -log10(obs.p), plot.it=F)

    # make a legend label from filename
    s <- unlist(strsplit(filename, "/"))
    lab <- unlist(strsplit(s[length(s)], "[.]"))[1]
    labels[i] <- lab

    if(i == 1) {
        plot(vals$x, vals$y, col=pal[i], las=1,
             xlab="null -log10 p-values",
             ylab="observed -log10 p-values")
        abline(a=0, b=1)
    } else {
        points(vals$x, vals$y, col=pal[i])
    }
}

legend("topleft", legend=labels, pch=20, col=pal[1:length(labels)])


dev.off()
