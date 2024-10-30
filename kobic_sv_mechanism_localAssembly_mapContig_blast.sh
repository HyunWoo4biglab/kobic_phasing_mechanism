#!/bin/bash


# Check if the input file argument is provided
if [ "$#" -eq 0 ]; then
    echo "Usage: $0 <inputdir> <working directory>"
    exit 1
fi

workdir="$1"
genome="$2"
bamdir="$3"

assembly_src="/home/hyunwoo/src/python/svProject/dnsv_workflow/sv_mechanism/dnsv_breakpoint_localAssembly_mapContig.sh"
extract_contig_src="/home/hyunwoo/src/python/svProject/dnsv_workflow/sv_mechanism/extract_velvet_contig_seq.py"

cd $workdir

# Check if the input directory exists
if [ ! -d "$workdir" ]; then
    echo "Error: Input directory not found."
    exit 1
fi

# Iterate over subdirectories under the input directory

for family in "$workdir"/*/; do
    # Check if the path is a directory
    if [ -d "$family" ]; then
        cd $workdir
        echo "Processing subdirectory: $family"
        targetdir="${family}/etching_call/"
        for subject in "$targetdir"/*/; do
            if [ -d "$subject" ]; then
                familynum=$(basename ${family})
                subjectnum=$(basename ${subject})
                if [[ "$subjectnum" != parent* ]]; then
                    echo "$family : $subject"
#                    target_regex=".+_${subjectnum}*.bam$"
                    sample_id=$(echo $subjectnum | awk -F "_" '{print $1}')
                    echo $sample_id
                    target_regex=".+_${sample_id}.final.bam$"
                    target_bam=$(find "$bamdir" -type l | grep -E "$target_regex")
                    #echo $target_bam
                    cd $subject
                    if [ ! -d "sv_mechanism" ]; then
                        mkdir "sv_mechanism"
                    fi
                    cd "sv_mechanism"
#                    outputdir="${subject}/sv_mechanism/"

                    targetfile="${subject}/${familynum}_${subjectnum}.scored.filtered.typed.contigFilt_50bp_PASS_nonLCR_nonBND.addVAF.dnsv.preciseModeFiltered.addHeader.callable_region_filtered.addClipRead.clipRead_filtered.vcf" ## check and update for dnSV vcf
                    while IFS= read -r line; do

                        match_ex="^#"
                        if [[ $line =~  $match_ex ]]; then
                            continue
                        fi

                        contig=$(echo "$line" | awk '{print $1}' | tr -d '\r')
                        start_pos=$(echo "$line" | awk '{print $2}' | tr -d '\r')
                        regex=";END=([0-9]+);"
                        if [[ "$line" =~ $regex ]]; then
                            end_pos="${BASH_REMATCH[1]}"
                        else
                            echo END POSITION NOT MATCHED!
                        fi
                        region="${contig}:${start_pos}-${end_pos}"

                        echo -e $subjectnum $region

                        if [ ! -d $region ]; then
                            mkdir $region
                        fi
                        cd $region

                        outputdir="${subject}/sv_mechanism/${region}/"
                        assembly_cmd="bash ${assembly_src} ${region} ${target_bam} ${outputdir} ${genome}"
                        echo $assembly_cmd
                        eval $assembly_cmd

                        #extract_output="${outputdir}/${region}_local_assembly_"
                        extract_cmd="python ${extract_contig_src} -i ${outputdir}/${region}_local_assembly_alignment.bam -r ${region} -o ${outputdir}/${region}_local_assembly_contig"
                        echo $extract_cmd
                        eval $extract_cmd

                        megablast_cmd="/home/hyunwoo/programs/blast/ncbi-blast-2.16.0+/bin/blastn -task megablast -query ${outputdir}/${region}_local_assembly_contig.bp1.fasta -subject ${outputdir}/${region}_local_assembly_contig.bp2.fasta -out ${outputdir}/${region}_local_assembly_contig.blastn_megablast_result.txt -outfmt 0"
                        echo $megablast_cmd
                        eval $megablast_cmd

                        smallblast_cmd="/home/hyunwoo/programs/blast/ncbi-blast-2.16.0+/bin/blastn -task blastn-short -query ${outputdir}/${region}_local_assembly_contig.bp1.fasta -subject ${outputdir}/${region}_local_assembly_contig.bp2.fasta -out ${outputdir}/${region}_local_assembly_contig.blastn_short_result.txt -outfmt 0"
                        echo $smallblast_cmd
                        eval $smallblast_cmd

                        cd ../

                    done < $targetfile
                fi
            fi
        done
    fi
done
