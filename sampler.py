import sys
import random
from Bio import SeqIO, motifs, Seq
from Bio.motifs import thresholds, matrix
import math
import random
import copy

nucleotides = ['A', 'C', 'G', 'T']
random.seed()

gapList={}


def main(filename, motiflength,gapsize):
    # read in file
    records = list(SeqIO.parse(filename, "fasta"))
    # creat random instances of given motif size
    instanceref = get_random_instances(records, motiflength,gapsize)
    print("found %d sequences" % len(instanceref))
    print("Got random instances:")
    motif = motifs.create(instanceref)
    print(motif)
    print("Starting profile")
    prof = create_pssm(motif)
    print(prof)
    motif = recursive_random(instanceref, motiflength, records,gapsize)
    print("Creating results.pdf")
    motif.weblogo("results.pdf", format="pdf", show_errorbars=False,
        show_ends=False, color_scheme="color_classic")


def create_pssm(instances):
    profile = instances.counts.normalize(pseudocounts=1)
    for nct in nucleotides:
        tup = []
        for i in range(len(profile[nct])):
            tup.append(-math.log(profile[nct][i]))
        profile[nct] = tuple(tup)
    return matrix.PositionSpecificScoringMatrix(alphabet=nucleotides,values=profile)


def recursive_random(instances, motiflength, records,gapsize):
    old_total = check_solution(instances)

    for idx, instance in enumerate(instances):
        train_instances = copy.deepcopy(instances)
        leave_out = random.choice(train_instances)
        seq_index=instances.index(leave_out)
        if (gapList[seq_index] == 0 or gapList[seq_index]==8):
            print("Leaving out %s" % leave_out)
        else:
            print("Leaving out %s" % leave_out[0:gapList[seq_index]] + "-" + leave_out[gapList[seq_index]:])
        train_instances.remove(leave_out)

        train_motifs = motifs.create(train_instances)
        profile = create_pssm(train_motifs)

        leftseqs = [records[seq_index]]
        new_instances = get_best_matches(leftseqs, profile, motiflength,gapsize,seq_index)
        print("new best instance:")
        for new_instance in new_instances:
            if(gapList[seq_index]==0 or gapList[seq_index]==8):
                print(new_instance)
            else:
                print((new_instance[0:gapList[seq_index]]+"-"+new_instance[gapList[seq_index]:]))
            instances[seq_index] = new_instance
    # printing the result from this iteration
    total = check_solution(instances)
    profile = create_pssm(motifs.create(instances))
    print("new solution: %d" % total)
    print("new profile: ")
    print(profile)
    #     Check if there is no regression, if not continue the recursion else stop the program
    if total < old_total:
        return recursive_random(instances,motiflength,records,gapsize)
    else:
        motif = motifs.create(instances)
        print("Finale Profile")
        print_pseudo(motif)
        print("Consensus sequence")
        print(motif.consensus)
        return motif


def print_pseudo(motif):
    counts = copy.deepcopy(motif.counts)
    for nct in nucleotides:
        tup = []
        for i in range(len(counts[nct])):
            tup.append(counts[nct][i]+1)
        counts[nct] = tuple(tup)
    print(counts)

def get_best_matches(leftseq, profile, motif_length,gapsize,idx):
    best_matches = list()
    for seq in leftseq:
        best_score = 999
        best_match = ""
        for i in range(len(seq) - 1 - motif_length-gapsize):
            for j in range(1,motif_length-2):
                # dnaseq = Seq.Seq(str(seq.seq)[i:i + motif_length]+)
                dnaseq = Seq.Seq(str(seq.seq)[i:i + j]+str(seq.seq)[i+j+gapsize:i+motif_length+gapsize])
                score = profile.calculate(dnaseq)
                if (score < best_score):
                    best_score = score
                    best_match = dnaseq
                    gapList[idx]=j
        best_matches.append(best_match)
    return best_matches


def check_solution(instances):
    pssm = create_pssm(motifs.create(instances))
    old_scores = list()
    for instance in instances:
        old_scores.append(pssm.calculate(instance))
    return sum(old_scores)


def get_random_instances(records, motiflength,gapsize):
    instances = motifs.Instances()
    for idx,record in enumerate(records):
        pos = random.randint(0, len(record.seq) - (motiflength+gapsize))
        gappos=random.randint(0,motiflength)
        seq=None
        if(pos+gappos+gapsize-pos+motiflength+gapsize>0):
            seq=record.seq[pos:pos + gappos]+record.seq[pos+gappos+gapsize:pos+motiflength+gapsize]
        else:
            if(gappos==0):
                seq=record.seq[pos + gappos + gapsize:pos + motiflength + gapsize]
            else:
                seq = record.seq[pos:pos + gappos]
        instances.append(seq)
        gapList[idx]=gappos
    return instances


if __name__ == "__main__":
    filename = ""
    motiflength = 0
    gapsize=0
    if len(sys.argv) < 2:
        print("usage: sampler.py [filename] [motiflength]")
    else:
        try:
            filename = sys.argv[1]
            motiflength = int(sys.argv[2])
            for idx,arg in enumerate(sys.argv[3:]):
                if(arg=="--gap"):
                    # + 4 as we start counting at idx 3 and need current index + 1
                    gapsize=int(sys.argv[idx+4])

        except:
            print("usage: sampler.py [filename] [motiflength]")
    main(filename, motiflength,gapsize)
