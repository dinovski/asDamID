# asDamID
allele-specific DamID sequencing pipeline

1. Trim raw reads for long/short damID adapters

2. Align trimmed reads with bowtie2 (default end-to-end, as required by SNPsplit) to a reference genome N-masked at high quality homozygous SNPs differing between parental strains (eg. CAST/129)

3. Assign allele-specific reads > spearate files for G1 specific, G2 specific, unassigned alignments, and conflictual.  
(Common alleles are those for which the maternal and paternal reads have the same number of mismatches and same position. If scores are the same but position is different between the reads, these ambiguous reads (UA) are output to a separate bam)  

    XX:Z:UA alignment is not assigned to any parental genome  
    XX:Z:G1 alignment is assigned to the first parental genome  
    XX:Z:G2 alignment is assigned to the second parental genome  
    XX:Z:CF alignment is ambiguous/conflictual

4. Filter reads that do not begin with the 'GATC' motif and MAPQ<10 > allele-specific alignment summary statistics

5. Generate bigWig files for log2 ratios and RPKM normalized counts

6. Divide reference genome coordinates into desired bin size (eg. 10kb)

7. Extract 'GATC' coordinates from the reference FASTA

8. Count reads mapping within 2 GATC sites and collapse into the genomic bin with the largest overlap

9. Normalize bins with > 10 counts to 1 million reads and calculated the log2 ratio of the fusion over dam counts using a pseudcount of 1. Compute the average ratio for replicates

10. Define LADs by running a Hidden Markov Model (https://github.com/gui11aume/HMMt) over the normalized counts
