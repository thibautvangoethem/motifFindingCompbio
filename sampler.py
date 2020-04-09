import sys
import random
from Bio import SeqIO, motifs
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
    for instance in instanceref:
        print(instance)
    profile = create_profile(instanceref)
    print("Starting profile")
    for nucleotide in nucleotides:
        print(nucleotide + "\t\t", end='')
        for instance in profile:
            # this special print syntax is to allign the printing
            print("{:>12} \t".format(str(profile[instance][nucleotide])), end='')
        # print a newline
        print("")
    recursive_random(instanceref, motiflength, records)


def create_profile(instances):
    pseudocount = 1
    instance_dict = {}
    for i in range(len(instances[0])):
        instance_dict[i] = {}
        instance_dict[i]["all"] = 0
        for nucleotide in nucleotides:
            instance_dict[i][nucleotide] = pseudocount
            instance_dict[i]["all"] += 1
        for seq in instances:
            instance_dict[i][seq[i]] += 1
            instance_dict[i]["all"] += 1
        for nucleotide in nucleotides:
            # We make a weighted calculation to get the profile numbers
            instance_dict[i][nucleotide] = -math.log(instance_dict[i][nucleotide] / instance_dict[i]["all"]);
            instance_dict[i][nucleotide] = int(instance_dict[i][nucleotide] * 100) / 100
    return instance_dict


def recursive_random(instances, motiflength, sequences):
    old_total = check_solution(instances)

    for idx, instance in enumerate(instances):
        leave_out = random.randint(0, len(instances) - 1)
        print("Leaving out %s" % instances[leave_out])
        # making training instances which are a copy of the given instances -1 randomly chosen instance
        train_instances = list()
        for train_idx, train_instance in enumerate(instances):
            if (train_idx != leave_out): train_instances.append(copy.deepcopy(train_instance))

        profile = create_profile(train_instances)

        leftseqs = [sequences[leave_out]]

        new_instances = get_best_matches(leftseqs, profile, motiflength)
        print("new best instance:")
        for new_instance in new_instances:
            print(new_instance)
            instances[leave_out] = new_instance
    # printing the result from this iteration
    total = check_solution(instances)
    profile = create_profile(instances)
    print("new solution: %d" % total)
    print("new profile: ")
    for nucleotide in nucleotides:
        print(nucleotide + "\t\t", end='')
        for instance in profile:
            # this special print syntax is to allign the printing
            print("{:>12} \t".format(str(profile[instance][nucleotide])), end='')
        # print a newline
        print("")
    #     Check if there is no regression, if not continue the recursion else stop the program
    if total < old_total:
        profile=recursive_random(instances,motiflength,sequences)
    else:
        pass
        # Here something strange is done in the perl script (with the scalar(@instances) and i do not know what it means
        """print("Final profile")
        for nucleotide in nucleotides:
            print(nucleotide + "\t\t", end='')
            for i in range(motiflength):
                # this special print syntax is to allign the printing
                print(str(int(math.exp(-profile[i][nucleotide]))))
                print("{:>12} \t".format(str(int(math.exp(-profile[i][nucleotide]))*((instances)+4))), end='')
            # print a newline
            print("")"""

def get_best_matches(leftseq, profile, motif_length):
    best_matches = list()
    for seq in leftseq:
        best_score = 999
        best_match = ""
        for i in range(len(seq) - 1 - motif_length):
            dnaseq = str(seq.seq)[i:i + motif_length]
            score = score_profile(dnaseq, profile)
            if (score < best_score):
                best_score = score
                best_match = dnaseq
        best_matches.append(best_match)
    return best_matches


def check_solution(instances):
    profile = create_profile(instances)
    old_scores = list()
    for instance in instances:
        old_scores.append(score_profile(instance, profile))

    old_total = 0
    for score in old_scores:
        old_total += score
    return old_total


def score_profile(seq, profile):
    score = 0
    for idx, nucleotide in enumerate(seq):
        score += profile[idx][nucleotide]
    return score


def get_random_instances(records, motiflength):
    instances = []
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
