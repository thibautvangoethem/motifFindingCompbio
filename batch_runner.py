import sampler
from config import Config
from jellyfish import jaro_distance
import sys
import os

import matplotlib.pyplot as plt

# This file is to use for testing/creating graphs only. Dont use it as an actuall sampler because you have less
# cotrol through parameters here
batchlist = [1, 4]
if __name__ == "__main__":
    filename = "data/complex_motif_more.fsa"
    motiflength = 8

    # type the motif that you know it will be here
    expected_motif = 'TATTAACA'
    #     everything else will be default  values from the config
    Distandvalue = []
    Dist = []
    for i in batchlist:
        Config.batch = i
        best = sampler.main(filename, motiflength)

        jarodist = jaro_distance(str(best[0].consensus), expected_motif)

        Dist.append(jarodist)
        Distandvalue.append([jarodist, best[0].consensus])
        sys.stdout = sys.__stdout__
        print(best[0].consensus)
        print(jarodist)
        sys.stdout = open(os.devnull, 'w')

    # printing the result
    # enable loggin for this part
    sys.stdout = sys.__stdout__
    print(Distandvalue)

    plt.plot(batchlist, Dist, 'ro')
    plt.show()
