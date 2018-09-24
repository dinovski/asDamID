#!/usr/bin/sh

## $1 = input file(s)
## $2 = output path
fastqc_func()
{
    local out=$2/fastqc
    local log=$3/fastqc.log
    mkdir -p ${out}
    
    local cmd="${FASTQC_PATH}/fastqc -t 4 -j ${JAVA_PATH}/java -o ${out} $1 2> $log"
    exec_cmd ${cmd} > ${out} 2>&1
}
