#!/bin/bash


# Check if the input file argument is provided
if [ "$#" -eq 0 ]; then
    echo "Usage: $0 <inputdir> <outputfile>"
    exit 1
fi

inputfile="$1"
workdir="$2"
bamdir="$3"
cohort="$4"
src="dnSV_parent_of_origin_predictor_240826.py" ####### change src code path
cd $workdir

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

    if [ $cohort == "CEPH" ]; then
#        child_regex="SRR[0-9]+_${child_id}_[^Child].+\.BQSR.Apply.bam$"
        child_regex="SRR[0-9]+_${child_id}_[cC]hild.+\.BQSR.Apply.bam$"
        paternal_regex="SRR[0-9]+_${paternal_id}_[^Child].+\.BQSR.Apply.bam$"
        maternal_regex="SRR[0-9]+_${maternal_id}_[^Child].+\.BQSR.Apply.bam$"
    elif [ $cohort == "ABC" ] || [ $cohort == 'CTC' ]; then
        child_regex=".+_${child_id}.final.bam$"
        paternal_regex=".+_${paternal_id}.final.bam$"
        maternal_regex=".+_${maternal_id}.final.bam$"
        #child_bai_regex=".+_${child_id}.final.bam.bai$"
        #paternal_bai_regex=".+_${paternal_id}.final.bam.bai$"
        #maternal_bai_regex=".+_${maternal_id}.final.bam.bai$"
    elif [ $cohort == "RD" ]; then
        child_regex=".+_${child_id}*.bam$"
        paternal_regex=".+_${paternal_id}*.bam$"
        maternal_regex=".+_${maternal_id}*.bam$"
    elif [ $cohort == "Radar" ]; then
        child_regex="${child_id}.*.bam$"
        paternal_regex="${paternal_id}.*.bam$"
        maternal_regex="${maternal_id}.*.bam$"
    fi

    child_bam=$(find "$bamdir" -type l | grep -E "$child_regex")
    paternal_bam=$(find "$bamdir" -type l | grep -E "$paternal_regex")
    maternal_bam=$(find "$bamdir" -type l | grep -E "$maternal_regex")

    cd $workdir


    family="family${family_id}"
    targetdir="${workdir}/${family}/etching_call/${child_id}/"

    targetfile="${targetdir}/${family}_${child_id}.scored.filtered.typed.contigFilt_50bp_PASS_nonLCR_nonIntraBND.addVAF.dnsv.nonBND.parentVCFfilter.clipRead_filter.vcf" ### check and change target file name
    output="${targetdir}/${family}_${child_id}_dnSV_phasing"

    if [ -f $targetfile ]; then
        cmd="python ${src} -v ${targetfile} -t sv -cb ${child_bam} -fb ${paternal_bam} -mb ${maternal_bam} -c 5 -o ${output}"
        echo $cmd
#        eval $cmd
    else
        echo "FILE NOT FOUND : ${targetfile}"
    fi

done < "$inputfile"
