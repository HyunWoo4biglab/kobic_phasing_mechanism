#!/bin/bash


# Check if the input file argument is provided
if [ "$#" -eq 0 ]; then
    echo "Usage: $0 <inputdir> <outputfile>"
    exit 1
fi

inputfile="$1"
workdir="$2"
bamdir="$3"

src="/home/hyunwoo/src/python/svProject/dnsv_workflow/sv_phasing/convert_vcf_to_unfazedBed.py" #change path

while IFS= read -r line; do
    ((line_num++))

    if [ $line_num -eq 1 ]; then
        continue
    fi

    family_id=$(echo "$line" | awk '{print $1}' | tr -d '\r')
    child_id=$(echo "$line" | awk '{print $2}' | tr -d '\r')
    paternal_id=$(echo "$line" | awk '{print $3}' | tr -d '\r')
    maternal_id=$(echo "$line" | awk '{print $4}' | tr -d '\r')

    child_regex="${child_id}*"
    paternal_regex="${paternal_id}*"
    maternal_regex="${maternal_id}*"

    child=$(find "$bamdir" -type l | grep -E "$child_regex")
    child_bam=$(echo $child | awk -F " " '{print $1}')
    if [[ $child_bam == *bai ]]; then
        child_bam=${child_bam%.bai}
        #child_bam=$(echo $child | awk -F " " '{print $2}')
    fi

    paternal=$(find "$bamdir" -type l | grep -E "$paternal_regex")
    paternal_bam=$(echo $paternal | awk -F " " '{print $1}')
    if [[ $paternal_bam == *bai ]]; then
        #echo $paternal_bam
        paternal_bam=${paternal_bam%.bai}
        #paternal_bam=$(echo $paternal | awk -F " " '{print $2}')
    fi

    maternal=$(find "$bamdir" -type l | grep -E "$maternal_regex")
    maternal_bam=$(echo $maternal | awk -F " " '{print $1}')
    if [[ $maternal_bam == *bai ]]; then
        maternal_bam=${maternal_bam%.bai}
        #maternal_bam=$(echo $maternal | awk -F " " '{print $2}')
    fi

    cd $workdir


    family="family${family_id}"
#    family="${family_id}"
    targetdir="${workdir}/${family}/etching_call/${family}_sample_${child_id}/"
#    targetdir="${workdir}/${family}/etching_call/${child_id}_maternal_${maternal_id}_paternal_${paternal_id}"

    cd $target_dir
    if [ ! -d phasing_unfazed ]; then
        mkdir phasing_unfazed
    fi

    targetfile="${targetdir}/${family}_sample_${child_id}.scored.filtered.typed.contigFilt_50bp_PASS_nonLCR_nonIntraBND.addVAF.dnsv.nonBND.parentVCFfilter.addClipRead.clipRead_filtered.vcf" ### check and change target file name
#    targetfile="${targetdir}/${family}_${child_id}_maternal_${maternal_id}_paternal_${paternal_id}.scored.filtered.typed.contigFilt_50bp_PASS_nonLCR_nonBND.addVAF.dnsv.preciseModeFiltered.addHeader.callable_region_filtered.addClipRead.clipRead_filtered.vcf"

    bed="${targetdir}/${family}_sample_${child_id}.dnsv_phasing_input.bed"
    ped="${targetdir}/${family}_sample_${child_id}.dnsv_phasing_input.ped"


    child_snv="child_snv" # NEED TO ADD PATH
    paternal_snv="paternal_snv" # NEED TO ADD PATH
    maternal_snv="maternal_snv" # NEED TO ADD PATH
    merged_snv="${targetdir}/${family}_${child_id}_${paternal_id}_${maternal_id}.merge.vcf.gz"
    unfazed_output="${targetdir}/${family}_${child_id}.unfazed_phasing.txt"

    if [ -f $targetfile ]; then
        # CONVERT DNSV VCF TO BED
        convert_to_bed_cmd="/BiO/scratch/users/minhak/anaconda3/envs/evaf/bin/python ${src} -v ${targetfile} -i ${child_id} -o ${bed}"
#        convert_to_bed_cmd="python ${src} -v ${targetfile} -i ${child_id} -o ${bed}"
        echo $convert_to_bed_cmd
#       eval $convert_to_bed_cmd

        # GENERATE FAMILY-MERGED SNV VCF
        compress_child="bgzip ${child_snv}"
        echo $compress_child
#        eval $compress_child

        compress_paternal="bgzip ${paternal_snv}"
        echo $compress_paternal
#        eval $compress_paternal

        compress_maternal="bgzip ${maternal_snv}"
        echo $compress_maternal
#        eval $compress_maternal

        index_child="tabix -p vcf ${child_snv}.gz"
        echo $index_child
#        eval $index_child

        index_paternal="tabix -p vcf ${paternal_snv}.gz"
        echo $index_paternal
#        eval $index_paternal

        index_maternal="tabix -p vcf ${maternal_snv}.gz"
        echo $index_maternal
#        eval $index_maternal

        make_merged_snv="bcftools merge -Oz ${child_snv}.gz ${paternal_snv}.gz ${maternal_snv}.gz -o ${merged_snv}"
        echo $make_merged_snv
#        eval $make_merged_snv
        
        make_ped_file="echo -e '${family}\t${child_id}\t${paternal_id}\t${maternal_id}\t0\t0' > ${ped}"
        echo $make_ped_file
        eval $make_ped_file

        # RUN UNFAZED
        unfazed_cmd="unfzed -d ${bed} -s ${merged_snv} -p ${ped} -g 38 --bam-pairs ${child_id}:${child_bam} ${paternal_id}:${paternal_bam} ${maternal_id}:${maternal_bam} > ${unfazed_output}"
        echo $unfazed_cmd
#        eval $unfazed_cmd

    else
        echo "FILE NOT FOUND : ${targetfile}"
    fi

done < $inputfile

