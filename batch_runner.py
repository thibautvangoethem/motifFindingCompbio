from sampler import *

# This file is to use for testing/creating graphs only. Dont use it as an actuall sampler because you have less cotrol through parameters here
batchlist=[1,4]
if __name__=="__main__":
    filename="data/complex_motif_more.fsa"
    motiflength=8
    # type the motif that you know it will be here
    expected_motif='TATTAACA'
#     everything else will be default  values from the config
    Distandvalue=[]
    Dist=[]
    for i in batchlist:
        sys.stdout = sys.__stdout__
        print(background(filename))
        motif_list = list()
        sys.stdout = open(os.devnull, 'w')

        for batch_counter in range(i):
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
            motifextra = recursive_random(instanceref, motiflength, records)  # Motifextra is a tuple of Motif, score and gaplength

            motif_list.append(motifextra)
        # getting the best motif
        best = None

        for motif in motif_list:
            if (best == None or motif[1] < best[1]): #The lower the score the better
                best = motif


        jarodist=jaro_distance(str(best[0].consensus), expected_motif)

        Dist.append(jarodist)
        Distandvalue.append([jarodist,best[0].consensus])
        sys.stdout = sys.__stdout__
        print(best[0].consensus)
        print(jarodist)
        sys.stdout = open(os.devnull, 'w')

    # printing the result
    # enable loggin for this part
    sys.stdout = sys.__stdout__
    print(Distandvalue)

    plt.plot(batchlist,Dist,'ro')
    plt.show()

    # plt.plot(Scoreval)
    # plt.show()

