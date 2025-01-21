#!/bin/bash


# Check if the input file argument is provided
if [ "$#" -eq 0 ]; then
    echo "Usage: $0 <input_pedigree> <vcf_directory> <output_directory>"
    exit 1
fi

inputfile="$1"
vcfdir="$2"
outputdir="$3"
genome="$4"
joint_outputdir="$5"
#cd $workdir

line_num=0
while IFS= read -r line; do
    ((line_num++))

    if [ $line_num -eq 1 ]; then
        continue
    fi

    family_id=$(echo "$line" | awk '{print $1}' | tr -d '\r')
    child_id=$(echo "$line" | awk '{print $2}' | tr -d '\r')
    paternal_id=$(echo "$line" | awk '{print $3}' | tr -d '\r')
    maternal_id=$(echo "$line" | awk '{print $4}' | tr -d '\r')
    contig=$(echo "$line" | awk '{print $5}' | tr -d '\r')
    start_pos=$(echo "$line" | awk '{print $6}' | tr -d '\r')
    end_pos=$(echo "$line" | awk '{print $7}' | tr -d '\r')
#    trio_type=$(echo "$line" | awk '{print $7}' | tr -d '\r')

    interval_start=$(($start_pos-100000))
    interval_end=$(($end_pos+100000))
    target_interval="${contig}:${interval_start}-${interval_end}"

    child_regex=".+_${child_id}*.g.vcf.gz$"
    paternal_regex=".+_${paternal_id}*.g.vcf.gz$"
    maternal_regex=".+_${maternal_id}*.g.vcf.gz$"

    child_vcf=$(find "$vcfdir" -type l | grep -E "$child_regex")
    #child_vcf="${vcfdir}/${child_id}/${child_id}.genotype.g.vcf.gz"
    paternal_vcf=$(find "$vcfdir" -type l | grep -E "$paternal_regex")
    #paternal_vcf="${vcfdir}/${paternal_id}/${paternal_id}}.genotype.g.vcf.gz"
    maternal_vcf=$(find "$vcfdir" -type l | grep -E "$maternal_regex")
    #maternal_vcf="${vcfdir}/${maternal_id}/${maternal_id}}.genotype.g.vcf.gz" 
    #/BiO/scratch/users/minhak/selee/output/HaplotypeCaller/concat_output/

    family="family${family_id}"
    trio_gdb_dir="${outputdir}/${family}_${child_id}_gdb/"
    tmp_dir="${trio_gdb_dir}/tmp/"
    if [ ! -d $trio_gdb_dir ]; then
        mkdir -p $tmp_dir
    fi

#    cd $trio_gdb_dir
#    if [ ! -d tmp ]; then
#        mkdir tmp
#    fi

    gdb_output="${trio_gdb_dir}/${family}_${child_id}_${paternal_id}_${maternal_id}_gdb_${contig}_${start_pos}"
    joint_output="${joint_outputdir}/${family}_${child_id}_${paternal_id}_${maternal_id}_${target_interval}.genotypes.vcf.gz"

    if [ -f $joint_output ]; then

        # BUILD GEOMICS DB
        gdb_cmd="gatk --java-options \"-Djava.io.tmpdir=${tmp_dir} -Xms2G -Xmx4G -XX:ParallelGCThreads=2\" GenomicsDBImport --genomicsdb-workspace-path ${gdb_output} -R ${genome} -V ${child_vcf} -V ${paternal_vcf} -V ${maternal_vcf} --tmp-dir ${tmp_dir} --intervals ${target_interval}"
        echo $gdb_cmd
#        eval $gdb_cmd

        # JOINT GENOTYPING
        joint_cmd="gatk --java-options \"-Djava.io.tmpdir=${tmp_dir} -Xms2G -Xmx2G -XX:ParallelGCThreads=2\" GenotypeGVCFs -R ${genome} -V gendb://${gdb_output} -O ${joint_output}"
        
    #    outputfile="${outputdir}/${family}_${child_id}_${paternal_id}_${maternal_id}.trioCombine.g.vcf"
    #    concat_cmd="gatk --java-options \"-Xmx4g\" CombineGVCFs -R ${genome} --variant ${child_vcf} --variant ${paternal_vcf} --variant ${maternal_vcf} -O ${outputfile}"
        echo $joint_cmd
#        eval $joint_cmd
    fi

    # UNFAZED
#    family="family${family_id}"
##    family="${family_id}"
#
#    targetdir="${workdir}/${family}/etching_call/${family}_sample_${child_id}/"
#
#    cd $target_dir
#    if [ ! -d phasing_unfazed ]; then
#        mkdir phasing_unfazed
#    fi
#
#    if [ $trio_type == "ABC_complete" ]; then    

done < $inputfile
