#
# Script to simulate data to be used with combined haplotype test
#

library(optparse)

option.list <- list(
  make_option(c("-a", "--alpha"),
              help="expected number of reads from REF chromosome overlapping a single site such as a heterozygous site",
              default=5.0, type="double"),
  
  make_option(c("-e", "--effect.size"), 
              help="ratio of reads from REF/NON-REF chromosome (alpha/beta) for 'true' sites", 
              default=1.2, type="double"),
  
  make_option(c("-r", "--region.size"), 
              help="regionsize * alpha is expected number of mapped reads in target region", 
              default=20, type="double"),
  
  make_option(c("-n", "--n.sample"),
              help="number of samples to generate data for",
              default=20, type="integer"),
  
  make_option(c("-m", "--n.site"),
              help="number of sites/regions to simulate data for",
              default=10000, type="integer"),
  
  make_option(c("-t", "--true.fraction"),
              help="fraction of 'true' sites to simulate that have non-zero effect size",
              default=0.1, type="double"),
  
  make_option(c("-f", "--allele.freq"),
              help="REF allele frequency in population", 
              default=0.1, type="double"),
  
  make_option(c("-o", "--output.dir"),
              help="directory to write output files to",
              default=".", type="character"),
  
  make_option(c("p", "--output.prefix"),
              help="prefix of output file",
              default="sim_hap_read_counts", type="character")
  
)

opt.parser <- OptionParser(option_list=option.list)

opt <- parse_args(opt.parser)

# flag whether site is from ALTERNATIVE (alpha != beta) or NULL (alpha == beta)
n.alt <- opt$n.site * opt$true.fraction
n.null <- opt$n.site * (1.0 - opt$true.fraction)
is.alt <-  c(rep(1, n.alt), rep(0, n.null))

site.names <- c(paste("ALT.", 1:n.alt, sep=""),
                paste("NULL.", 1:n.null, sep=""))

# beta depends on whether site is from ALT or NULL model 
# NULL sites have no difference between alleles (alpha==beta)
alpha <- rep(opt$alpha, n.alt + n.null)
beta <- alpha / ((opt$effect.size * is.alt) + (1.0-is.alt))

output.prefix <- paste(opt$output.dir, "/", opt$output.prefix, sep="")

# simulate genotypes at each site, given minor allele freq
# and assuming HWE
ref.af <- opt$allele.freq
expect.geno.freq <- c(ref.af*ref.af, 2*ref.af*(1-ref.af), (1-ref.af) * (1-ref.af))

# draw samples for each site, giving number of homozygous ref, het, homozygous alt individuals
genos <- rmultinom(opt$n.site, size=opt$n.sample, prob=expect.geno.freq)

# create a genotype matrix with samples as columns, rows as sites, using sample to shuffle order of genotypes
geno.matrix <- t(apply(genos, 2, function(x) { sample(c(rep(0, x[1]), rep(1, x[2]), rep(2, x[3])))}))

expect.ref.reads <- alpha * opt$region.size
expect.alt.reads <- beta * opt$region.size
expect.homo.ref.reads <- 2 * expect.ref.reads
expect.het.reads <- expect.ref.reads + expect.alt.reads
expect.homo.alt.reads <- 2 * expect.alt.reads

expect.het.site.reads <- alpha + beta
expect.het.site.ref.prp <- alpha / (alpha + beta)

# loop over each sample... generating read counts and allele-specific read counts
for(col in 1:ncol(geno.matrix)) {
   gtypes <- geno.matrix[,col]
   
   homo.ref.idx <- which(gtypes == 0)
   het.idx <- which(gtypes == 1)
   homo.alt.idx <- which(gtypes == 2)
   
   n.homo.ref <- length(homo.ref.idx)
   n.het <- length(het.idx)
   n.homo.alt <- length(homo.alt.idx)
   
   # generate counts from each allele... 
   # for read counts this should be scaled up for region size
   # TODO: replace poisson with beta-negative binomial
   homo.ref.reads <- rpois(n.homo.ref, expect.homo.ref.reads)
   het.reads <- rpois(n.het, expect.het.reads)
   homo.alt.reads <- rpois(n.homo.alt, expect.homo.alt.reads)
   read.counts <- rep(NA, opt$n.site)
   read.counts[homo.ref.idx] <- homo.ref.reads
   read.counts[het.idx] <-het.reads
   read.counts[homo.alt.idx] <- homo.alt.reads
   
   het.site.reads <- rpois(n.het, expect.het.site.reads)
   
   
   # generate counts overlapping heterozygous site
   # currently assume single heterzygous site...
   # TODO: replace binom with beta-binom
   ref.allele.count <- rep(NA, length(gtypes))
   alt.allele.count <- rep(NA, length(gtypes))

   ref.allele.count[het.idx] <- rbinom(n.het, het.site.reads, expect.het.site.ref.prp)
   alt.allele.count[het.idx] <- het.site.reads - ref.allele.count[het.idx]
   
   hap.str <- rep("", opt$n.site)
   hap.str[homo.ref.idx] <- "0|0"
   hap.str[het.idx] <- "0|1"
   hap.str[homo.alt.idx] <- "1|1"
   
   het.prob <- rep(0.01, opt$n.site)
   het.prob[het.idx] <- 0.99
   
   # create an output table for each sample
   out.tab <- data.frame(CHROM=rep("chr1", opt$n.site),
                        TEST.SNP.POS=rep(1000, opt$n.site),
                        TEST.SNP.ID=site.names,
                        TEST.SNP.REF.ALLELE=rep("A", opt$n.site),
                        TEST.SNP.ALT.ALLELE=rep("G", opt$n.site),
                        TEST.SNP.GENOTYPE=gtypes,
                        TEST.SNP.HAPLOTYPE=hap.str,
                        REGION.START=rep(1, opt$n.site),
                        REGION.END=rep(2000, opt$n.site),
                        REGION.SNP.POS=rep(1000, opt$n.site),
                        REGION.SNP.HET.PROB=het.prob,
                        REGION.SNP.LINKAGE.PROB=rep(1.0, opt$n.site),
                        REGION.SNP.REF.HAP.COUNT=ref.allele.count,
                        REGION.SNP.ALT.HAP.COUNT=alt.allele.count,
                        REGION.SNP.OTHER.HAP.COUNT=rep(0, opt$n.site),
                        REGION.READ.COUNT=read.counts,
                        GENOMEWIDE.READ.COUNT=rep(10e6, opt$n.site))
              
   filename <- paste(output.prefix, "_", col, ".txt", sep="")
   write.table(out.tab, file=filename, quote=F, sep="\t", row.names=F, col.names=T)
}



