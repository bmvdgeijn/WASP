
base.dir <- "/iblm/netapp/home/gmcvicker/data1/gmcvicker/data/WASP_power_calcs"

power.sim.tab <- read.table(paste(base.dir, "/power_sim_params.txt", sep=""), header=T)

tab <- power.sim.tab
tab$N.SIGNIF <- rep(NA, nrow(tab))
tab$N.TEST <- rep(NA, nrow(tab))
tab$FRAC.SIGNIF <- rep(NA, nrow(tab))

for(i in seq(1, nrow(tab))) {
  ref.af <- tab$REF.ALLELE.FREQ[i]
  effect.size <- tab$EFFECT.SIZE[i]
  samp.size <- tab$SAMPLE.SIZE[i]
  
  data.dir <- sprintf("%s/sim_af%.2f_es%.2f_ss%d", base.dir, ref.af, effect.size, samp.size)

  # results.as.file <- paste(data.dir, "/cht_results.as_only.txt", sep="")
  # results.bnb.file <- paste(data.dir, "/cht_results.bnb_only.txt", sep="")
  results.combined.file <- paste(data.dir, "/cht_results.txt", sep="")

  if(file.exists(results.combined.file)) {
    # as.tab <- read.table(results.as.file, header=T)
    # bnb.tab <- read.table(results.bnb.file, header=T)
    combined.tab <- read.table(results.combined.file, header=T)

    n.test <- nrow(combined.tab)
    expect.p <- seq(1, n.test)/(1 + n.test)
    obs.p <- combined.tab$P.VALUE
    #qqplot(-log10(expect.p), -log10(obs.p))
    #abline(a=0, b=1, col="red")

    # use Benjamini-Hochberg method to calculate FDR
    fdr <- p.adjust(obs.p, method="BH")

    # what p-value corresponds to FDR 5%?
    if(any(fdr < 0.05)) {
      p.thresh <- max(obs.p[fdr <= 0.05])
      # lines(x=c(0, log10(n.test)), y=rep(-log10(p.thresh), 2), col="blue", lty=2)

      n.signif <- sum(obs.p <= p.thresh)
      frac.signif <- n.signif / n.test
    } else {
      n.signif <- 0
      frac.signif <- 0.0
    }
  
      tab$N.SIGNIF[i] <- n.signif
      tab$N.TEST[i] <- n.test
      tab$FRAC.SIGNIF[i] <- frac.signif
  }
}

library(RColorBrewer)


effect.sizes <- unique(tab$EFFECT.SIZE)
samp.sizes <- seq(10, 60, by=10)
col.pal <- brewer.pal(length(samp.sizes), "Set1")

pdf("cht_power_analysis.pdf", width=5, height=10)

par(mfrow=(c(2,1)))

for(j in 1:length(samp.sizes)) {
  f <- tab$EFFECT.SIZE == 1.2 & 
      tab$SAMPLE.SIZE == samp.sizes[j] & 
      tab$REF.ALLELE.FREQ <= 0.5
  
  if(j == 1) {
      plot(tab$REF.ALLELE.FREQ[f], tab$FRAC.SIGNIF[f], col=col.pal[j], ylim=c(0,1),
           type="o", las=1, ylab="prob signif test at FDR < 0.05", pch=20, 
           main=paste("log2(effect size)=0.26", sep=""), xlab="minor allele freq")
      
      legend("bottomright", legend=paste("n=", samp.sizes, sep=""), 
             pch=20, col=col.pal, title="sample size")
    } else {
      points(tab$REF.ALLELE.FREQ[f], tab$FRAC.SIGNIF[f], pch=20, col=col.pal[j], type="o")
    }
  }


for(i in 1:length(samp.sizes)) {
  f <- tab$SAMPLE.SIZE == samp.sizes[i] & tab$REF.ALLELE.FREQ == 0.20 & abs(log2(tab$EFFECT.SIZE)) <= 1
  
  new.tab <- tab[f,]
  new.tab <- new.tab[order(new.tab$EFFECT.SIZE),]
  
  if(i == 1) {
    plot(log2(new.tab$EFFECT.SIZE), new.tab$FRAC.SIGNIF, col=col.pal[i], ylim=c(0,1),
         type="o", las=1, ylab="prob signif test at FDR < 0.05", pch=20, 
         main=paste("minor allele freq=", 0.2, sep=""), xlab="log2(effect size)")
    
    legend("bottomright", legend=paste("n=", samp.sizes, sep=""), 
           pch=20, col=col.pal, title="sample size")
  } else {
    points(log2(new.tab$EFFECT.SIZE), new.tab$FRAC.SIGNIF, pch=20, col=col.pal[i], type="o")
  }
}

dev.off()





