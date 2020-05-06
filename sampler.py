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
gapSize=0
maxGapsize=0


def main(filename, motiflength):
    # read in file
    records = list(SeqIO.parse(filename, "fasta"))
    # creat random instances of given motif size
    instanceref = get_random_instances(records, motiflength)
    print("found %d sequences" % len(instanceref))
    print("Got random instances:")
    motif = motifs.create(instanceref)
    print(motif)
    print("Starting profile")
    prof = create_pssm(motif)
    print(prof)
    motif = recursive_random(instanceref, motiflength, records)
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


def recursive_random(instances, motiflength, records):
    global gapSize
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
        new_instances = get_best_matches(leftseqs, profile, motiflength,seq_index)
        print("new best instance:")
        for new_instance in new_instances:
            if(gapList[seq_index]==0 or gapList[seq_index]==8):
                print(new_instance)
            else:
                print("gapsize: "+str(gapSize))
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
        return recursive_random(instances,motiflength,records)
    else:
        motif = motifs.create(instances)
        print("Finale Profile")
        print_pseudo(motif)
        print("Consensus sequence")
        print("gapsize: " + str(gapSize))
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

def get_best_matches(leftseq, profile, motif_length,idx):
    global gapSize
    global gapList
    best_matches_scores=list()
    inner_Gaplist=copy.deepcopy(gapList)
    # We loop over all the possible gapsizes to search for the one with the best total_score
    for loop_gapsize in range(maxGapsize+1):
        total_score=0
        best_matches = list()
        # Loop over all sequences and get the best match
        for seq in leftseq:
            best_score = 999
            best_match = ""
            # Loop over all motifs in the sequence
            for i in range(len(seq) - 1 - motif_length-loop_gapsize):
                for j in range(1,motif_length-2):
                    # This long line is for getting a sequence with a given gap in it
                    dnaseq = Seq.Seq(str(seq.seq)[i:i + j]+str(seq.seq)[i+j+loop_gapsize:i+motif_length+loop_gapsize])
                    score = profile.calculate(dnaseq)
                    if (score < best_score):
                        best_score = score
                        inner_Gaplist[idx]=j
                        best_match = (dnaseq)
            best_matches.append(best_match)
            total_score+=best_score
        best_matches_scores.append((total_score,best_matches,copy.deepcopy(inner_Gaplist)))
    best= best_matches_scores[0]
    # get the best gap size sequence
    for idx,match in enumerate(best_matches_scores):
        if(match[0]<best[0]):
            best=match
            gapSize=idx
    gapList=best[2]
    return best[1]


def check_solution(instances):
    pssm = create_pssm(motifs.create(instances))
    old_scores = list()
    for instance in instances:
        old_scores.append(pssm.calculate(instance))
    return sum(old_scores)

def get_random_instances(records, motiflength):
    # get a random gapsize to start, the gap will be refined after multiple iterations of the algorithm
    global gapSize
    gapSize=random.randint(0,maxGapsize+1)
    instances = motifs.Instances()
    for idx,record in enumerate(records):
        pos = random.randint(0, len(record.seq) - (motiflength+gapSize))
        gappos=random.randint(0,motiflength)
        seq=None
        if(pos+gappos+gapSize-pos+motiflength+gapSize>0):
            seq=record.seq[pos:pos + gappos]+record.seq[pos+gappos+gapSize:pos+motiflength+gapSize]
        else:
            if(gappos==0):
                seq=record.seq[pos + gappos + gapSize:pos + motiflength + gapSize]
            else:
                seq = record.seq[pos:pos + gappos]
        instances.append(seq)
        gapList[idx]=gappos
    return instances


if __name__ == "__main__":
    filename = ""
    motiflength = 0
    if len(sys.argv) < 2:
        print("usage: sampler.py [filename] [motiflength]")
    else:
        try:
            filename = sys.argv[1]
            motiflength = int(sys.argv[2])
            for idx,arg in enumerate(sys.argv[3:]):
                if(arg=="--maxgap"):
                    # + 4 as we start counting at idx 3 and need current index + 1
                    maxGapsize=int(sys.argv[idx+4])

        except:
            print("usage: sampler.py [filename] [motiflength]")
    main(filename, motiflength)
