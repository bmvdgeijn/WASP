import genome.db, gzip, argparse, math

def main():
    error=0.01
    args=parse_options()
    if args.infile[-3:]==".gz":
        infile=gzip.open(args.infile,"r")
    else:
        infile=open(args.infile,"r")
    if args.outfile[-3:]==".gz":
        outfile=gzip.open(args.outfile,"w")
    else:
        outfile=open(args.outfile,"w")

    gdb = genome.db.GenomeDB()
    ref_track=gdb.open_track(args.ref_track)
    alt_track=gdb.open_track(args.alt_track)

    snp_line=infile.readline()
    if snp_line:
        outfile.write(snp_line)
    else:
        sys.stderr.write("The input file was empty.\n")
        exit()

    snp_line=infile.readline()
    while snp_line:
        snpinfo=snp_line.strip().split()
        if snpinfo[9]=="NA":
            outfile.write(snp_line)
        else:
            new_hetps=process_one_snp(snpinfo, ref_track, alt_track,error)
            outfile.write("\t".join(snpinfo[:10]+[";".join(new_hetps)]+snpinfo[11:])+"\n")
        snp_line=infile.readline()
    

def process_one_snp(snpinfo, ref_track, alt_track,error):
    chrm=snpinfo[0]
    # positions of target SNPs
    snplocs=[int(y.strip()) for y in snpinfo[9].split(';')]
    # heterozygote probabilities of target SNPs
    hetps = [float(y.strip()) for y in snpinfo[10].split(';')]
    update_hetps=[]
    for i in range(len(snplocs)):
        pos=snplocs[i]
        adr=ref_track.get_nparray(chrm, pos, pos)[0]
        ada=alt_track.get_nparray(chrm, pos, pos)[0]
        update_hetps.append(str(get_posterior_hetp(hetps[i],adr,ada,error)))
    return update_hetps


def get_posterior_hetp(hetp_prior,adr,ada,error):
    prior = min(.99,hetp_prior)
    badlike=addlogs(math.log(error)*adr+math.log(1-error)*ada,math.log(1-error)*adr+math.log(error)*ada)
    goodlike=math.log(0.5)*adr+math.log(0.5)*ada
    if goodlike-badlike > 40:
        return 1
    else:
        return prior*math.exp(goodlike-badlike)/(prior*math.exp(goodlike-badlike)+(1-prior))
        
def addlogs(loga,logb):
    return max(loga,logb)+math.log(1+math.exp(-abs(loga-logb)))

def parse_options():
    parser=argparse.ArgumentParser()
    parser.add_argument("infile", action='store', default=None)
    parser.add_argument("outfile", action='store', default=None)
    parser.add_argument("ref_track", action='store', default=None)
    parser.add_argument("alt_track", action='store', default=None)
    
    return parser.parse_args()

main()
