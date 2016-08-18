import sys
import string
import subprocess
import os


DNA_COMP = None

def comp(seq_str):
    """complements the provided DNA sequence and returns it"""
    global DNA_COMP

    if DNA_COMP is None:
        DNA_COMP = string.maketrans("ATCGMRWSYKNatcgmrwsykn",
                                    "TAGCKYWSRMNtagckywsrmn")
    return seq_str.translate(DNA_COMP)


def revcomp(seq_str):
    """returns reverse complement of provided DNA sequence"""
    return comp(seq_str)[::-1]

        
def sort_bam(input_bam, output_prefix):
    """Calls samtools sort on input_bam filename and writes to
    output_bam. Takes into account that the command line arguments 
    for samtools sort have changed between versions."""

    output_bam = output_prefix + ".sort.bam"
    
    # first try new way of using samtools sort
    failed = False
    cmd = "samtools sort -o " + output_bam + " " + input_bam
    sys.stderr.write("running command: %s\n" % cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except Exception as e:
        sys.stderr.write("samtools sort command failed:\n%s\n" %
                         str(e))
        failed = True
    if not os.path.exists(output_bam):
        sys.stderr.write("output file %s does not exist\n" % output_bam)
        failed = True
        
    if failed:
        # OLD way of calling samtools (changed in newer versions)
        sys.stderr.write("samtools sort command failed, trying old samtools "
                         "syntax\n")
        
        cmd = "samtools sort " + input_bam + " " + output_prefix
        sys.stderr.write("running command: %s\n" % cmd)

        try:
            subprocess.check_call(cmd, shell=True)
        except Exception as e:
            sys.stderr.write("samtools sort command failed:\n%s\n" %
                             str(e))
            exit(1)
        
        if not os.path.exists(paths.sorted_output_bam):
            raise IOError("Failed to create sorted BAM file '%s'" %
                          paths.sorted_output_bam)
