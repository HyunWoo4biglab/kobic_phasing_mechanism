#!/bin/bash


# Check if the input file argument is provided
if [ "$#" -eq 0 ]; then
    echo "Usage: $0 <inputdir> <outputfile>"
    exit 1
fi

inputfile="$1"
#workdir="$2"
#bamdir="$3"
vcfdir="$2"
#cohort="$5"


while IFS= read -r line; do
    ((line_num++))

    if [ $line_num -eq 1 ]; then
        continue
    fi

    family_id=$(echo "$line" | awk '{print $1}' | tr -d '\r')
    child_id=$(echo "$line" | awk '{print $2}' | tr -d '\r')
    paternal_id=$(echo "$line" | awk '{print $4}' | tr -d '\r')
    maternal_id=$(echo "$line" | awk '{print $3}' | tr -d '\r')
    contig=$(echo "$line" | awk '{print $5}' | tr -d '\r')
    start_pos=$(echo "$line" | awk '{print $6}' | tr -d '\r')
    end_pos=$(echo "$line" | awk '{print $7}' | tr -d '\r')
    class=$(echo "$line" | awk '{print $8}' | tr -d '\r')

    if [ $class == "complete" ]; then
        workdir="/home/hyunwoo/data_set/project/svProject/wgs/abc/complete_trio/"
        bamdir="/home/hyunwoo/data_set/project/svProject/wgs/abc/bam_symbolic_abc_total/"
#        vcfdir="/home/hyunwoo/data_set/project/svProject/wgs/abc/gvcf_abc_total/"

    elif [ $class == "reverse_father" ]; then
        workdir="/home/hyunwoo/data_set/project/svProject/wgs/abc/incomplete_trio/reverse_trio_father/sv_call/"
        bamdir="/home/hyunwoo/data_set/project/svProject/wgs/abc/bam_symbolic_abc_total/"
#        vcfdir="/home/hyunwoo/data_set/project/svProject/wgs/abc/gvcf_abc_total/"


    elif [ $class == "reverse_mother" ]; then
        workdir="/home/hyunwoo/data_set/project/svProject/wgs/abc/incomplete_trio/reverse_trio_mother/sv_call/"
        bamdir="/home/hyunwoo/data_set/project/svProject/wgs/abc/bam_symbolic_abc_total/"
#        vcfdir="/home/hyunwoo/data_set/project/svProject/wgs/abc/gvcf_abc_total/"
    fi

    child_regex=".+_${child_id}.final.bam$"
    paternal_regex=".+_${paternal_id}.final.bam$"
    maternal_regex=".+_${maternal_id}.final.bam$"
    re_start=$((start_pos-100000))
    re_end=$((end_pos+100000))
    joint_vcf="${vcfdir}/family${family_id}_${child_id}_${paternal_id}_${maternal_id}_${contig}:${re_start}-${re_end}.genotypes.vcf.gz"

#    echo $joint_regrex

#    joint_vcf=$(find "$vcfdir" -type f | grep -e "$joint_regrex")

    bed="${workdir}/Family${family_id}/etching_call/${child_id}_maternal_${maternal_id}_paternal_${paternal_id}/phasing_unfazed/Family${family_id}_${child_id}.dnsv_phasing_input.bed"
    ped="${workdir}/Family${family_id}/etching_call/${child_id}_maternal_${maternal_id}_paternal_${paternal_id}/phasing_unfazed/Family${family_id}_${child_id}.dnsv_phasing_input.ped"


    child_bam=$(find "$bamdir" -type l | grep -E "$child_regex")
    #child_bam=$(echo $child | awk -F " " '{print $1}')
    #if [[ $child_bam == *bai ]]; then
    #    child_bam=${child_bam%.bai}
    #    #child_bam=$(echo $child | awk -F " " '{print $2}')
    #fi

    paternal_bam=$(find "$bamdir" -type l | grep -E "$paternal_regex")
    #paternal_bam=$(echo $paternal | awk -F " " '{print $1}')
    #if [[ $paternal_bam == *bai ]]; then
    #    #echo $paternal_bam
    #    paternal_bam=${paternal_bam%.bai}
    #    #paternal_bam=$(echo $paternal | awk -F " " '{print $2}')
    #fi

    maternal_bam=$(find "$bamdir" -type l | grep -E "$maternal_regex")
    #maternal_bam=$(echo $maternal | awk -F " " '{print $1}')
    #if [[ $maternal_bam == *bai ]]; then
    #    maternal_bam=${maternal_bam%.bai}
    #    #maternal_bam=$(echo $maternal | awk -F " " '{print $2}')
    #fi


    child_basename=$(basename $child_bam)
    paternal_basename=$(basename $paternal_bam)
    maternal_basename=$(basename $maternal_bam)
    
    child_prefix="${child_basename%.final.bam}"
    maternal_prefix="${maternal_basename%.final.bam}"
    paternal_prefix="${paternal_basename%.final.bam}"

    outputfile="${workdir}/Family${family_id}/etching_call/${child_id}_maternal_${maternal_id}_paternal_${paternal_id}/phasing_unfazed/Family${family_id}_${child_id}.unfazed_phasing.txt"
    cmd="unfazed -d ${bed} -s ${joint_vcf} -p ${ped} -g 38 --bam-pairs ${child_prefix}:${child_bam} ${paternal_prefix}:${paternal_bam} ${maternal_prefix}:${maternal_bam} > ${outputfile}"
    echo $cmd
    eval $cmd

done < $inputfile

