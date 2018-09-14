#!/usr/bin/env python

from Bio import SeqIO
import re
import os
import sys
import csv

usage="""
Perform case insensitive motif search of a fasta file and output 2 bed files: motif_sites.bed (motif coordinates) and motif_fragments.bed (coordinates between each motif)
python fasta2bed.py <motif> <fasta.fa> <out_path>
"""

if len(sys.argv) < 4:
    print usage

pattern = sys.argv[1]
fasta_path = sys.argv[2]
out_path = sys.argv[3]
motif_sites = os.path.join(out_path, pattern + "_sites.bed")
motif_frags = os.path.join(out_path, pattern + "_fragments.bed")

outfile_sites = open(motif_sites, 'w')

def search_fasta(pattern, fasta_path):
    for record in SeqIO.parse(open(fasta_path, 'rU'), "fasta"):
       chrom = record.id
       for match in re.finditer(pattern, str(record.seq), re.IGNORECASE):
           start_pos = match.start() + 1
           end_pos = match.end()
           outfile_sites.write(chrom + "\t" + str(start_pos) + "\t" + str(end_pos) + "\n")

search_fasta(pattern, fasta_path)
outfile_sites.close()

print "converting sites to fragments"

## convert sites to fragments: coordinates are from motif1 start to motif2 end
infile = motif_sites
outfile_frags = open(motif_frags, 'w')

with open(infile,'r') as f:
    reader = csv.reader(f, delimiter='\t')
    previous = reader.next()
    
    for line in reader:
        if line[0] == previous[0]: #to process by chr
            chr = previous[0]
            start = int(previous[1])
            end = int(line[2])
            #print chr, start, end
            previous = line
            outfile_frags.write(chr + "\t" + str(start) + "\t" + str(end) + "\n")
        else:
            print "processing %s" % (line[0])
            previous = line

outfile_frags.close()

if __name__ == "__main__":
    main()
