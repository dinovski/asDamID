#!/bin/bash

CUTADAPT=/bioinfo/local/build/Centos/python/python-2.7.12/bin/cutadapt
FASTX_REVCOM=/bioinfo/local/build/fastx_toolkit_0.0.13/fastx_reverse_complement
TRIMG=/bioinfo/local/build/Centos/trim_galore/trim_galore_v0.4.4

IDIR=/data/kdi_prod/.kdi/project_workspace_0/1242/acl/01.00/svn/analyse/script/damID
ODIR=/data/tmp/dina/damid/trimmed
FILE=/data/kdi_prod/.kdi/project_workspace_0/1242/acl/01.00/svn/analyse/script/damID/FASTQS_NEW
FASTQ=$(grep "^$PBS_ARRAYID," $FILE | cut -d',' -f2)

mkdir -p ${ODIR}

InPutBaseName=`basename ${FASTQ}`

# Set some parameters
MIN_Q=25

# set base filename
case $InPutBaseName in
  *.fq.gz)
    OutFileName=${InPutBaseName%.fq.gz}
    CAT=zcat ;;
  *.fastq.gz)
    OutFileName=${InPutBaseName%.fastq.gz}
    CAT=zcat ;;
  *.fastq)
    OutFileName=${InPutBaseName%.fastq}
    CAT=cat ;;
  *.fq)
    OutFileName=${InPutBaseName%.fq}
    CAT=cat ;;
  *)
  echo "inputfile (${InPutBaseName}) should either be gzipped fastq"
  echo "(*.fq.gz or fastq.gz) or fastq (*.fq or *.fastq)"
  exit 1 ;;
esac

#####
## Adaptor sequences used in DamID-seq
#####
ADPTR_SHORT_5="GGTCGCGGCCGAG"
ADPTR_LONG_5="CTAATACGACTCACTATAGGGCAGCGTGGTCGCGGCCGAG"

# reverse complement of adapter sequences
#  (${FASTX_REVCOM} expects fasta input, "awk 'NR > 1'"
ADPTR_SHORT_3=`echo -e ">\n${ADPTR_SHORT_5}" | ${FASTX_REVCOM} | awk 'NR > 1'`
ADPTR_LONG_3=`echo -e ">\n${ADPTR_LONG_5}" | ${FASTX_REVCOM} | awk 'NR > 1'`

CLIP_STATS=${OutFileName}_local.clip_stats

###############################################################################
## trim adapter sequences from reads, split fastq in inner and edge reads
###############################################################################
## clip long adapter, save non-adapter reads, and adapter-match stats
ADPTR_LONG_LEN=`echo "${ADPTR_LONG_5}" | wc -c`
ADPTR_LONG_LEN=`expr ${ADPTR_LONG_LEN} - 1`
LEN_THRES=`expr ${ADPTR_LONG_LEN} / 2`

# Cutadapt
#  -g. 5' adapter
#  -a. 3' adapter
#  -O. Minimum overlap length between the read and the adapter (in order to trim the overlap)
#  --match-read-wildcards. Enables to consider IUPAC degenerated nucleotides in the reads
#  --discard-trimmed. Throw away reads in which an adapter was found. (In this case
#  -o. output File

# Discard reads that have at least 20bp of the adapters
OutLongUnTrimmed=${ODIR}/untrimmed-long_${OutFileName}.fq
OutLongTrimmed=${ODIR}/trimmed-long_${OutFileName}.fq
$CUTADAPT -g $ADPTR_LONG_5 -a $ADPTR_LONG_3 -O $LEN_THRES -m 20 --match-read-wildcards ${FASTQ} --untrimmed-output ${OutLongUnTrimmed} -o ${OutLongTrimmed} > ${ODIR}/$CLIP_STATS

# Clip short adapter, save non-adapter and trimmed reads in separate files, and save adapter-match stats
OutShortUnTrimmed=${ODIR}/untrimmed-short_${OutFileName}.fq
OutShortTrimmed=${ODIR}/trimmed-short_${OutFileName}.fq
cat ${OutLongUnTrimmed} | $CUTADAPT -g $ADPTR_SHORT_5 -a $ADPTR_SHORT_3 -m 20 --match-read-wildcards - --untrimmed-output $OutShortUnTrimmed -o $OutShortTrimmed >> ${ODIR}/$CLIP_STATS

#cat ${TMP_FQ_INNER} | (${BOWTIE2} ${BOWTIE_PAR} -U - | samtools view -bS - -o ${TMP_BAM_INNER} ) 2> ${TMP_STATS_INNER}
#cat ${TMP_FQ_EDGE} | (${BOWTIE2} ${BOWTIE_PAR} -U - | samtools view -bS - -o ${TMP_BAM_EDGE} ) 2> ${TMP_STATS_EDGE}

## combine all long trimmed and short trimmed/untrimmed
cat ${OutLongTrimmed} ${OutShortTrimmed} ${OutShortUnTrimmed} | gzip -c > ${ODIR}/${OutFileName}_trimmed.fastq.gz
rm -rf ${OutLongTrimmed} ${OutLongUnTrimmed} ${OutShortTrimmed} ${OutShortUnTrimmed}

exit 0
