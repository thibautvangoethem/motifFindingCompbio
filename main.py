from Bio import SeqIO

#Class and method documentation: http://biopython.org/DIST/docs/api/Bio-module.html
#User friendly documentation: https://biopython.org/wiki/Documentation

for record in SeqIO.parse("data/complex_motif_more.fsa", "fasta"):
    print(record.id)
    print(record.seq)
