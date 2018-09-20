# asDamID
allele-specific DamID sequencing pipeline

1. Trim raw reads using cutadapt

2. Align trimmed reads with bowtie2 (default end-to-end, as required by SNPsplit) to mm10 N-masked at high quality homozygous SNPs differing between parental strains (eg. CAST/129)

3. Assign allele-specific reads using SNPsplit and split into parental allele specific BAM files
Common alleles are those for which the maternal and paternal reads have the same number of mismatches and same position. If scores are the same but position is different between the reads, these ambiguous reads (UA) are output to a separate bam.  

    XX:Z:UA alignment is not assigned to any parental genome  
    XX:Z:G1 alignment is assigned to the first parental genome  
    XX:Z:G2 alignment is assigned to the second parental genome  
    XX:Z:CF alignment is ambiguous/conflictual

4. Generate bigWig files for log2 ratios and RPKM normalized counts

5. Divide reference genome coordinates into desired bin size (eg. 10kb) using BEDTools windowMaker 

6. Extract 'GATC' coordinates from the reference FASTA

7. Reads mapping within 2 GATC sites are counted using featureCounts --largestOverlap and collapsed into the bin with the largest overlap (intersecBed -f 0.51 with binned genome)

8. Normalize bins with > 10 counts to 1 million reads and calculated the log2 ratio of the fusion over dam counts using a pseudcount of 1. Compute the average ratio for replicates

9. Define LADs by running a Hidden Markov Model (https://github.com/gui11aume/HMMt) over the normalized counts
