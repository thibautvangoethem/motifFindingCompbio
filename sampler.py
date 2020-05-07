
import sys

from Bio import SeqIO, motifs, Seq
from Bio.motifs import thresholds, matrix
from jellyfish import jaro_distance
import math
import random
import copy


from Bio.SeqRecord import SeqRecord

nucleotides = ['A', 'C', 'G', 'T']
random.seed()
palindromic = ["-p", "--palindrome"]

gapList={}
gapSize=0
maxGapsize=0


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
    Backgroundvector=background("Controledata.fsa") #Creert een vector met de relatieve frequentie van elke letter

    dumj=0 # Ik weet niet hoe ik anders over Backgroundvector moet lopen :/
    for nct in nucleotides:
        tup = []
        for i in range(len(profile[nct])):
            # Deel door het relevante gewicht van backgroundvector, als er meer A's zijn delen we door een getal groter dan 1
            #Hierdoor daalt het gewicht van een A.
            tup.append(-math.log((profile[nct][i])/Backgroundvector[dumj]))
        profile[nct] = tuple(tup)
        dumj=dumj+1
    return matrix.PositionSpecificScoringMatrix(alphabet=nucleotides,values=profile)


def recursive_random(instances, motiflength, records, palindrome=False):
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
        return recursive_random(instances, motiflength, records, palindrome)
    else:
        motif = motifs.create(instances)
        print("Finale Profile")
        print_pseudo(motif)
        print("Consensus sequence")
        print("gapsize: " + str(gapSize))
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

# def get_best_matches(leftseq, profile, motif_length,idx):
#     global gapSize
#     global gapList
#     best_matches_scores=list()
#     inner_Gaplist=copy.deepcopy(gapList)
#     # We loop over all the possible gapsizes to search for the one with the best total_score
#     for loop_gapsize in range(maxGapsize+1):
#         total_score=0
#         best_matches = list()
#         # Loop over all sequences and get the best match
#         for seq in leftseq:
#             best_score = 999
#             best_match = ""
#             # Loop over all motifs in the sequence
#             for i in range(len(seq) - 1 - motif_length-loop_gapsize):
#                 for j in range(1,motif_length-2):
#                     # This long line is for getting a sequence with a given gap in it
#                     dnaseq = Seq.Seq(str(seq.seq)[i:i + j]+str(seq.seq)[i+j+loop_gapsize:i+motif_length+loop_gapsize])
#                     score = profile.calculate(dnaseq)
#                     if (score < best_score):
#                         best_score = score
#                         inner_Gaplist[idx]=j
#                         best_match = (dnaseq)
#             best_matches.append(best_match)
#             total_score+=best_score
#         best_matches_scores.append((total_score,best_matches,copy.deepcopy(inner_Gaplist)))
#     best= best_matches_scores[0]
#     # get the best gap size sequence
#     for idx,match in enumerate(best_matches_scores):
#         if(match[0]<best[0]):
#             best=match
#             gapSize=idx
#     gapList=best[2]
#     return best[1]


def check_solution(instances):
    pssm = create_pssm(motifs.create(instances))
    old_scores = list()
    for instance in instances:
        old_scores.append(pssm.calculate(instance))
    return sum(old_scores)


def get_random_instances(records, motiflength, palindrome=None):
    # get a random gapsize to start, the gap will be refined after multiple iterations of the algorithm
    global gapSize
    gapSize=random.randint(0,maxGapsize+1)
    instances = motifs.Instances()
    for idx, record in enumerate(records):
        pos = random.randint(0, len(record.seq) - (motiflength + gapSize))
        gappos = random.randint(0, motiflength)
        seq = None
        if (pos + gappos + gapSize - pos + motiflength + gapSize > 0):
            seq = record.seq[pos:pos + gappos] + record.seq[pos + gappos + gapSize:pos + motiflength + gapSize]
        else:
            if (gappos == 0):
                seq = record.seq[pos + gappos + gapSize:pos + motiflength + gapSize]
            else:
                seq = record.seq[pos:pos + gappos]
        instances.append(seq)
        gapList[idx] = gappos
    return instances

#Een functie die zelf een motief genereert in een achtergrond
def create_control_data(Stringlength=70,Motif='TATTAA',amountofstrings=8,Errorpercentage=0):
    Recordlist = amountofstrings * [0] #Een lijst met seqrecord objecten

    Motif = list(Motif) #Makkelijker te manipuleren dan strings


    for i in range(len(Recordlist)):
        Background = round(Stringlength / 4) * ['A', 'C', 'G', 'T'] #Homogene verdeling
        random.shuffle(Background) #Randomiseer
        Background = Background[0:Stringlength] #Zorgt ervoor dat niet alle strings 4*n lang moeten zijn
        beginpoint = random.randint(0, Stringlength - len(Motif)) #Willekeurige locatie voor het begin van het motief
        Motif_mutated = Motif
        if random.random() < Errorpercentage / 100 * 1.25:  # Kan muteren naar dezelfde letter daarom 25 % meer kans dan aangegeven

            Motif_mutated[random.randint(0, len(Motif) - 1)] = random.choice(nucleotides)  # Kan muteren naar dezelfde letter

        Promotor = Background[0:beginpoint] + Motif_mutated + Background[beginpoint + len(Motif):] #Het hele stuk DNA

        Promotorstr = ''.join(Promotor) #Als string ipv list voor seq.seq functie

        Recordlist[i] = SeqRecord(Seq.Seq(Promotorstr), id='Fake' + str(i), description='Part of test data ') #Creert een SeqRecord object

    with open("Controledata.fsa", "w") as output_handle:
        SeqIO.write(Recordlist, output_handle, "fasta") #Schrijft alles naar een FASTA File
    return Recordlist


def background(filename):
    records = list(SeqIO.parse(filename, "fasta"))
    Totalsequence = Seq.Seq('')

    for recseq in records:
        Sequence = recseq.seq
        Totalsequence += Sequence #Neemt alle sequenties samen

    Backgroundfrequency = []
    for letter in nucleotides:
        Backgroundfrequency.append(Totalsequence.count(letter) / len(Totalsequence)) #Berekent de frequentie


    factor_vector = []
    for Freq in Backgroundfrequency:
        factor_vector.append(Freq / min(Backgroundfrequency)) #Berekent de relatieve frequentie

    return factor_vector

def usage():
    print("usage: sampler.py [option] [filename] [motiflength]")
    print("Options:")
    print("\t -p [float] --palindrome [float]\n\t\tTry to find palindromic motif with ratio higher as [float]")
    print("\t\tRatio is defined as the jaro distance between a sequence and its reverse complement")
    exit(0)

if __name__ == "__main__":
    Controldata = create_control_data(150, 'TATTAACCA', 15,0)
    pal = False
    # for i in range(len(Controldata)):
    #     print(Controldata[i].seq) #Geeft alle gebruikte strings weer

    filename = ""

    motiflength = 0

    if len(sys.argv) < 3:
        usage()
    else:
        try:
            if len(sys.argv) > 3:
                if sys.argv[1] in palindromic:
                    pal = float(sys.argv[2])
                    filename = sys.argv[3]
                    motiflength = int(sys.argv[4])
                    #
                    # filename = sys.argv[1]
                    # motiflength = int(sys.argv[2])
                    # for idx, arg in enumerate(sys.argv[3:]):
                    #     if (arg == "--maxgap"):
                    #         # + 4 as we start counting at idx 3 and need current index + 1
                    #         maxGapsize = int(sys.argv[idx + 4])
            else:
                filename = sys.argv[1]
                motiflength = int(sys.argv[2])
        except Exception as e:
            usage()
    main(filename, motiflength, palindrome=pal)




