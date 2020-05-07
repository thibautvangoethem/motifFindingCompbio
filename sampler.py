import sys
import random
from Bio import SeqIO, motifs, Seq
from Bio.motifs import thresholds, matrix
from jellyfish import jaro_distance
import math
import random
import copy

nucleotides = ['A', 'C', 'G', 'T']
random.seed()
palindromic = ["-p", "--palindrome"]


def main(filename, motiflength, palindrome=None):
    # read in file
    records = list(SeqIO.parse(filename, "fasta"))
    # creat random instances of given motif size
    instanceref = get_random_instances(records, motiflength, palindrome)
    print("found %d sequences" % len(instanceref))
    print("Got random instances:")
    motif = motifs.create(instanceref)
    print(motif)
    print("Starting profile")
    prof = create_pssm(motif)
    print(prof)
    motif = recursive_random(instanceref, motiflength, records, palindrome)
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


def recursive_random(instances, motiflength, records, palindrome=False):
    old_total = check_solution(instances)

    for idx, instance in enumerate(instances):
        train_instances = copy.deepcopy(instances)
        leave_out = random.choice(train_instances)
        print("Leaving out %s" % leave_out)
        train_instances.remove(leave_out)

        train_motifs = motifs.create(train_instances)
        profile = create_pssm(train_motifs)

        leftseqs = [records[instances.index(leave_out)]]
        new_instances = get_best_matches(leftseqs, profile, motiflength, palindrome)
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
        return recursive_random(instances, motiflength, records, palindrome)
    else:
        motif = motifs.create(instances)
        print("Finale Profile")
        print_pseudo(motif)
        print("Consensus sequence")
        print("new solution: %d" % total)
        if palindrome:
            print(motif.consensus + "----" + motif.consensus.reverse_complement())
            print(jaro_distance(str(motif.consensus), str(motif.consensus.reverse_complement())))
        else:
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


def get_best_matches(leftseq, profile, motif_length, palindrome):
    best_matches = list()
    for seq in leftseq:
        best_score = 999
        best_match = ""
        for i in range(len(seq) - 1 - motif_length):
            dnaseq = Seq.Seq(str(seq.seq)[i:i + motif_length])
            rev_compl = dnaseq.reverse_complement()
            if jaro_distance(str(dnaseq), str(rev_compl)) > palindrome:
                score = profile.calculate(dnaseq)
                if (score < best_score):
                    best_score = score
                    best_match = dnaseq
        if best_match == "":
            exit("No sequence found with a high enough palindromic ratio, try lowering it")
        best_matches.append(best_match)
    return best_matches


def check_solution(instances):
    pssm = create_pssm(motifs.create(instances))
    old_scores = list()
    for instance in instances:
        old_scores.append(pssm.calculate(instance))
    return sum(old_scores)


def get_random_instances(records, motiflength, palindrome=None):
    instances = motifs.Instances()
    for record in records:
        """if palindrome:
            found = False
            while not found:
                pos = random.randint(0, len(record.seq) - motiflength)
                seq = record.seq[pos:pos + motiflength]
                rev_compl = seq.reverse_complement()
                if jaro_distance(seq, rev_compl) > palindrome:
                    instances.append(seq)
                    found = True
        else:"""
        pos = random.randint(0, len(record.seq) - motiflength)
        seq = record.seq[pos:pos + motiflength]
        instances.append(seq)
    return instances

def usage():
    print("usage: sampler.py [option] [filename] [motiflength]")
    print("Options:")
    print("\t -p [float] --palindrome [float]\n\t\tTry to find palindromic motif with ratio higher as [float]")
    print("\t\tRatio is defined as the jaro distance between a sequence and its reverse complement")
    exit(0)

if __name__ == "__main__":
    filename = ""
    motiflength = 0
    pal = False
    if len(sys.argv) < 3:
        usage()
    else:
        try:
            if len(sys.argv) > 3:
                if sys.argv[1] in palindromic:
                    pal = float(sys.argv[2])
                    filename = sys.argv[3]
                    motiflength = int(sys.argv[4])
            else:
                filename = sys.argv[1]
                motiflength = int(sys.argv[2])
        except Exception as e:
            usage()
    main(filename, motiflength, palindrome=pal)
