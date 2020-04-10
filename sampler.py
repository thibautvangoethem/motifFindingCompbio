import sys
import random
from Bio import SeqIO, motifs, Seq
from Bio.motifs import thresholds, matrix
import math
import random
import copy

nucleotides = ['A', 'C', 'G', 'T']
random.seed()


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
    old_total = check_solution(instances)

    for idx, instance in enumerate(instances):
        train_instances = copy.deepcopy(instances)
        leave_out = random.choice(train_instances)
        print("Leaving out %s" % leave_out)
        train_instances.remove(leave_out)

        train_motifs = motifs.create(train_instances)
        profile = create_pssm(train_motifs)

        leftseqs = [records[instances.index(leave_out)]]
        new_instances = get_best_matches(leftseqs, profile, motiflength)
        print("new best instance:")
        for new_instance in new_instances:
            print(new_instance)
            instances[instances.index(leave_out)] = new_instance
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

def get_best_matches(leftseq, profile, motif_length):
    best_matches = list()
    for seq in leftseq:
        best_score = 999
        best_match = ""
        for i in range(len(seq) - 1 - motif_length):
            dnaseq = Seq.Seq(str(seq.seq)[i:i + motif_length])
            score = profile.calculate(dnaseq)
            if (score < best_score):
                best_score = score
                best_match = dnaseq
        best_matches.append(best_match)
    return best_matches


def check_solution(instances):
    pssm = create_pssm(motifs.create(instances))
    old_scores = list()
    for instance in instances:
        old_scores.append(pssm.calculate(instance))
    return sum(old_scores)


def get_random_instances(records, motiflength):
    instances = motifs.Instances()
    for record in records:
        pos = random.randint(0, len(record.seq) - motiflength)
        instances.append(record.seq[pos:pos + motiflength])
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
