
import sys
import numpy as np

# default values
def_allele_freq = 0.2
def_effect_size = 1.2
def_sample_size = 50

sys.stdout.write("ID REF.ALLELE.FREQ EFFECT.SIZE SAMPLE.SIZE\n")

sim_id = 0

                                           
# range of sample sizes
samp_size_step = 10
for samp_size in range(10, 60 + samp_size_step, samp_size_step):
    # sim_id += 1
    # sys.stdout.write("%d %.2f %.2f %d\n" %
    #                  (sim_id, def_allele_freq, def_effect_size, samp_size))

    # range of allele freqs
    step = 0.05
    for af in np.arange(0.05, 0.50 + step, step):
        sim_id += 1
        sys.stdout.write("%d %.2f %.2f %d\n" %
                         (sim_id, af, def_effect_size, samp_size))

    # range of effect sizes
    step = 0.2
    for log2_effect_size in np.arange(-1, 1 + step, step):
        effect_size = 2**log2_effect_size
        sim_id += 1
        sys.stdout.write("%d %.2f %.2f %d\n" %
                         (sim_id, def_allele_freq, effect_size, samp_size))
