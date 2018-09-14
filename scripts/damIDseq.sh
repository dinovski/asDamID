#!/bin/bash

## post alignment processing

## bin paths
SAMTOOLS=
BEDTOOLS=
PYTHON=
FEATURECOUNTS=
UCSC_TOOLS=
PYTHON=
DATAMASH=
R=

IDIR=
ODIR=
SCRIPTS_PATH=

LEN=mm10.chrom.sizes
FASTA=mm10.fa
BLACKLIST=mm10.blacklist.chrM.bed
BIN_SIZE=10000

## 01 prepare gatc fragment and binned genome bed/saf files
grep -v "chrUn" | grep -v "random" | grep -v "chrM" ${LEN} > ${IDIR}/chrom.sizes.filt
LEN=${IDIR}/chrom.sizes.filt

## get genome-wide GATC coordinates and fragments
${PYTHON} ${SCRIPTS_PATH}/fasta2bed.py gatc ${FASTA} ${IDIR}

## bed to SAF: only chr1-19, X, Y; remove sites that overlap (ie. GATC sites in close proximity)
grep -v "chrUn" ${ODIR}/gatc_fragments.bed | grep -v "_random" | grep -v "chrM" | sort -k1,1 -k2n,2 | \
awk '($3 > $2)' | awk -v OFS='\t' '{print "gatc"NR, $0, "+"}' > ${ODIR}/gatc_fragments.saf

## bin genome: only chr1-19, X, Y
${BEDTOOLS}/windowMaker -g ${LEN} -w ${BIN_SIZE} | awk -v OFS='\t' '{print "bin"NR, $0, "."}' > ${ODIR}/mm10_bins.saf

awk -v OFS='\t' '{print $2,$3,$4,$1}' ${ODIR}/mm10_bins.saf > ${ODIR}/mm10_bins.bed

## 02 generate binned counts and merge fusion/dam results
## count methods:
## gatc: only count reads mapping fully within a gatc fragment then collapse into bins
## bins: assign reads to bins with the largest overlap; still problem of gatc sites being in separate bins; need to look into this
COUNT_METHOD=gatc

GATC_SAF=${ODIR}/gatc_fragments.saf
BIN_BED=${ODIR}/mm10_bins.bed
BIN_SAF=${ODIR}/mm10_bins.saf

FC_OPTS='--largestOverlap'

## exclude blacklist regions
#BAM=${DIR}/${NAME}/mapping_N-masked_bowtie2/${NAME}_N-masked.allele_flagged.sort.bam
#${BEDTOOLS}/intersect -a ${BAM} -b ${BLACKLIST} -v > ${DIR}/${NAME}/mapping_N-masked_bowtie2/${NAME}.filt.bam

## count function: per gatc or bin, sort by binID
## output is bin_id,chr,start,end,count
bin_counts() {
	FILT_BAM=$1
	ODIR=$2
	if [ ${COUNT_METHOD} == "gatc" ];then
		${FEATURECOUNTS} ${FC_OPTS} -F SAF -a ${GATC_SAF} ${FILT_BAM} -o ${ODIR}/${NAME}_counts.txt
		## gatc counts per bin: require 51% overlap
		sed 1,2d ${ODIR}/${NAME}_counts.txt | awk -v OFS='\t' '{print $2,$3,$4,$7}' |\
		${BEDTOOLS}/intersectBed -f 0.51 -wo -a - -b ${BIN_BED} | awk -v OFS='\t' '{print $5,$6,$7,$4,$8}' |\
		${DATAMASH} -s -g 1,2,3,5 sum 4 | awk -v OFS='\t' '{print $4,$1,$2,$3,$5}' | sort -k1V,1 > ${ODIR}/${NAME}_counts.tab
	elif [ ${COUNT_METHOD} == "bins" ];then
		${FEATURECOUNTS} ${FC_OPTS} -F SAF -a ${BIN_SAF} ${FILT_BAM} -o ${ODIR}/${NAME}_counts.txt
		sed 1,2d ${ODIR}/${NAME}_counts.txt | awk -v OFS='\t' '{print $1,$2,$3,$4,$7}' | sort -k1V,1 > ${ODIR}/${NAME}_counts.tab
	else
		die "count method should be by 'gatc' or 'bins'"
	fi
}

echo "counting reads by $COUNT_METHOD"
echo "Total counts"
TOT_BAM=${DIR}/${NAME}/mapping_N-masked_bowtie2/${NAME}_N-masked.allele_flagged.sort.bam
TOT_ODIR=${DIR}/counts
mkdir -p ${TOT_ODIR}

bin_counts ${TOT_BAM} ${TOT_ODIR}

G1_BAM=${ODIR}/${NAME}/mapping_N-masked_bowtie2/G1/${NAME}_N-masked.genome1.sort.bam
G1_ODIR=${ODIR}/counts_g1
mkdir -p ${G1_ODIR}

bin_counts ${G1_BAM} ${G1_ODIR}

G2_BAM=${ODIR}/${NAME}/mapping_N-masked_bowtie2/G1/${NAME}_N-masked.genome2.sort.bam
G2_ODIR=${ODIR}/counts_g2
mkdir -p ${G2_ODIR}

bin_counts ${G2_BAM} ${G2_ODIR}

## 03 merge fusion/dam results

merge_counts () {
	ID=$1
	DAM=$2
	LMN=$3
	ODIR=$4
	for (( i=1; i<=$TOT; i++ )); do
		LC_ALL=C join -j 1 -t $'\t' ${ODIR}/${LMN}_counts.tab ${ODIR}/${DAM}_counts.tab |\
		awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$9}' > ${ODIR}/${ID}_merged_counts.txt
	done
}

G1_DIR=${DIR}/counts_g1
G2_DIR=${DIR}/counts_g2
merge_counts ${NAME} ${DAM} ${LMN} ${G1_DIR}
merge_counts ${NAME} ${DAM} ${LMN} ${G2_DIR}

## 04 combine replicates, sort by coordinate
LC_ALL=C join -j 1 -t $'\t' ${IDIR}/${NAME}_1_merged_counts.txt ${IDIR}/${NAME}_2_merged_counts.txt | sort -k2V,2 -k3n,3 |\
awk -v OFS="\t" 'BEGIN{OFS="\t"; print "BIN_ID","CHR","START","END","FUSION_1","DAM_1","FUSION_2","DAM_2"} ; {print $1,$2,$3,$4,$5,$6,$10,$11}' > ${IDIR}/${NAME}_merged_replicates.txt
done

## 05 run HMMt
cmd="${R}/R --no-save --no-restore -e \"scripts='${SCRIPTS_PATH}'; inPath='${IDIR}'; outPath='${ODIR}'; min_reads=10; source('${SCRIPTS_PATH}/hmmt.R', chdir=TRUE))\""
exec_cmd ${cmd} > $logdir/report.log 2>&1








