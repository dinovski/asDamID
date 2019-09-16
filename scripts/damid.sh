#!/usr/bin/sh

## $1 = output dir
## $2 = log dir
build_genome_func()
{
    local log=$2/build_genome.log
    local out=$1/genome
    mkdir -p ${out}

    local cmd="${SNPSPLIT}/SNPsplit_genome_preparation --nmasking --strain ${ID_REF} --strain2 ${ID_ALT} --reference_genome ${REF_FASTA} --vcf_file ${VCF}
    FASTA_NMASK=${out}/N-masked_${ID_REF}_${ID_ALT}.fa 2> $log"
    exec_cmd ${cmd} > ${out} 2>&1
}

## $1 = input file(s)
## $2 = output dir
## $3 = log dir
fastqc_func()
{
    local out=$2/fastqc
    local log=$3/fastqc.log
    mkdir -p ${out}
    
    local cmd="${FASTQC_PATH}/fastqc -t 4 -j ${JAVA_PATH}/java -o ${out} $1 2> $log"
    exec_cmd ${cmd} > ${out} 2>&1
}

nmasked_mapping_func() {

}

## $1 = input file(s)
## $2 = output dir
## $3 = log dir
markdup_func()
{
    local out=$2/mapping
    local log=$3/markdup.log
    out_bam=$(basename $1 ".bam")
    
    local cmd="${JAVA} -jar ${PICARD} MarkDuplicates \
    INPUT=$1 \
    OUTPUT=${out}/${out_bam}_noDup.bam \
    REMOVE_DUPLICATES=false \
    ASSUME_SORTED=true \
    CREATE_INDEX=true \
    METRICS_FILE=${out}/{out_bam}_metrics"
    
    exec_cmd ${cmd}
}

assign_reads_func()
{

}
