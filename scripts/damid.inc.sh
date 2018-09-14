## allele-specific DamID functions
## Copyright (c) 2018 Institut Curie
## Author(s): Dina Zielinski
## Contact: dina.zielinski@curie.fr
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.

## $1 = input file(s)
## $2 = output path
## $3 = log dir
fastqc_func()
{
    local out=$2/fastqc
    local log=$3/fastqc.log
    mkdir -p ${out}
    
    echo -e "Running FastQC ..."
    echo -e "Logs: ${log}"
    echo

    local cmd="${FASTQC_PATH}/fastqc -t 4 -j ${JAVA_PATH}/java -o ${out} $1 2> $log"
    exec_cmd ${cmd} > ${log} 2>&1
}
