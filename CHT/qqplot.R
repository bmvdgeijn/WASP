
# read output filename and  list of input filenames containing 
# CHT results from command line
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  usage <- "Usage:\n RScript --vanilla OUTPUT.png CHT_OUTPUT_1.txt [... CHT_OUTPUT_N.txt]"
  stop(paste("At least two arguments must be supplied.\n", usage), call.=FALSE)
}

png.filename <- args[1]
input.filenames <- args[2:length(args)]

# check that output filename looks like a PNG file
if(length(grep(".png$", png.filename)) == 0) {
  stop("expected output filename to end with .png", call.=FALSE)
}

# choose set of colors
library(RColorBrewer)
pal <- brewer.pal(9, "Set1")

png(png.filename, width=500, height=500)

labels <- c()

min.p <- 1e-20

for(i in 1:length(input.filenames)) {
    filename <- input.filenames[i]

    tab <- read.table(filename, header=T)
    n.test <- nrow(tab)
    null.p <- (1:n.test)/(n.test)
    obs.p <- tab$P.VALUE

    # cap p-values at min.p for drawing purposes
    obs.p[obs.p < min.p]  <- min.p
    null.p[null.p < min.p] <- min.p

    vals <- qqplot(-log10(null.p), -log10(obs.p), plot.it=F)

    # make a legend label from filename, stripping off extension and leading directories
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

