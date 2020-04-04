import sys
import random
from Bio import SeqIO, motifs

def main(filename, motiflength):
    records = SeqIO.parse("data/complex_motif_more.fsa", "fasta")
    instanceref = get_random_instances(records,motiflength)


def get_random_instances(records, motiflength):
    instances = []
    for record in records:
        pos = random.randint(0,len(record.seq)-motiflength)
        instances.append(record.seq[pos:pos+motiflength])
    return instances


if __name__ == "__main__":
    filename = ""
    motiflength = 0
    if len(sys.argv) != 3:
        print("usage: sampler.py [filename] [motiflength]")
    else:
        try:
            filename = sys.argv[1]
            motiflength = int(sys.argv[2])
        except:
            print("usage: sampler.py [filename] [motiflength]")
    main(filename, motiflength)
