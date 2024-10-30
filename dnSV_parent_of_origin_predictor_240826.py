#####################
# Written by mhchoi #
# Date: 2024-07-15  #
#####################

# REVISED BY HYUN WOO KIM
# 08-26-2024


import sys, argparse, subprocess, pysam, os, math, gzip 
from multiprocessing import Process, Queue

def variant_makeDic(vcf_file):
    check_vDic = dict(); var_count = 0

    if '.gz' not in vcf_file:
        f = open(vcf_file, 'r')
#        f = open(vcf_file, 'r')
    else:
        f = gzip.open(vcf_file, 'r')

    for lines in f:
        if lines.startswith('#'):   continue
        items = lines.strip().split('\t')
        var_count += 1

        contig = items[0]; start = int(items[1]); end = int(re.search(r'END=(\d+);', items[7]).group(1))
        svtype = re.search(r'SVTYPE=([A-Z]+);', items[7]).group(1)
        pos = f'{start}-{end}'

        # Format: {chrN: {Pos: {Ref: Variants}}}
        # Format: {chrN: {Pos: {Alt: list(Variants ..)}}}
        # format : {contig : {start_pos-end_pos : {sv_type : SVTYPE}}}
        if not contig in check_vDic:
            check_vDic.setdefault(contig, dict())#.setdefault(pos, dict()).setdefault('type', svtype)
            #check_vDic[contig][start].setdefault('end', list()).append(end)

        if not pos in check_vDic[contig]:
            check_vDic[contig].setdefault(pos, dict())
        if not 'type' in check_vDic[contig][pos]:
            check_vDic[contig][pos]['type'] = svtype

    f.close()

    print ('# This sample has {0} dnSVs.'.format(var_count))
    print(check_vDic)

    return check_vDic, var_count


def candidate_Vars_parser(read_info, pos_info, phasing_candidate_Dic):
    ref_positions = read_info.get_reference_positions()
    #ref_positions = read_info.get_reference_positions(full_length=True)
    ref_sequence = read_info.get_reference_sequence()
    q_sequence = read_info.query_sequence

    for q_pos in range(read_info.query_alignment_start, read_info.query_alignment_end):
        var_index = q_pos - read_info.query_alignment_start

        if var_index < len(ref_positions): # NEED EXPLANATION
#        if var_index < pos_info:
            #print(read_info.reference_name, pos_info)
            #print(var_index, len(ref_positions))
            ref_pos = ref_positions[var_index]

            # SKIP dnSV STARTING POSITION
            if int(ref_pos+1) != pos_info:
                if var_index < len(ref_sequence):
                    ref_base = ref_sequence[var_index]
                    q_base = q_sequence[q_pos]

                    # Ref와 다른 것을 확인함.
                    if ref_base.upper() != q_base.upper():
                        var_line = '{0}-{1}'.format(ref_base.upper(), q_base.upper())

                        # Save format: {Variant_position: {Ref-Alt: supporting read Count}}
                        if ref_pos not in phasing_candidate_Dic:
                            phasing_candidate_Dic.setdefault(ref_pos, dict()).setdefault(var_line, 1)
                        else:
                            if var_line not in phasing_candidate_Dic[ref_pos]:
                                phasing_candidate_Dic[ref_pos].setdefault(var_line, 1)
                            else:
                                phasing_candidate_Dic[ref_pos][var_line] += 1

                    else: pass
                else: pass
            else: pass
        else: pass

    return phasing_candidate_Dic

def fetch_mate(pos_info, bamfile, read_info, phasing_candidate_Dic):
    mateSam_file = pysam.AlignmentFile(bamfile, 'rb', threads=3)

    # WHY SUBSTRACTING 3?
    for mp_line in mateSam_file.pileup(read_info.next_reference_name, \
                                        read_info.next_reference_start-3, \
                                        read_info.next_reference_start+1, \
                                        stepper='all', min_mapping_quality=20, \
                                        ignore_overlaps=False, ignore_orphans=False, \
                                        min_base_quality=20):

        for pileupread in mp_line.pileups:
            # 오직 SNV type만 고려한다면?
            if (not pileupread.is_del) and (not pileupread.is_refskip):
                mate_info = pileupread.alignment

                if (mate_info.query_name == read_info.query_name) and (mate_info.is_read1 != read_info.is_read1):
                    phasing_candidate_Dic = candidate_Vars_parser(mate_info, pos_info, phasing_candidate_Dic)
                else: pass
            else: pass

    mateSam_file.close()

    return phasing_candidate_Dic

#################### checked



def check_only_DNV_reads(phasing_candidate_Dic, non_phasing_Dic):
    rm_Vars_Dic = dict()

    # DNV read만 가지는 특징이어야 하지 않을까?
    for tmp_pos, tmp_Var_Dic in phasing_candidate_Dic.items():
        for tmp_Var_line in tmp_Var_Dic.keys():
            if (tmp_pos in non_phasing_Dic) and (tmp_Var_line in non_phasing_Dic[tmp_pos]):
                rm_Vars_Dic.setdefault(tmp_pos, set()).add(tmp_Var_line)
            else: pass

    for tmp_pos, tmp_Var_set in rm_Vars_Dic.items():
        del phasing_candidate_Dic[tmp_pos]    

    return phasing_candidate_Dic

#################### checked 240827

def collect_the_variants_in_specific_position(samfile, chr_name, can_ref_pos, can_Var_list):
    tmp_candidate_dic = dict()
    # WHY can_ref_pos-3
    for mp_line in samfile.pileup(chr_name, can_ref_pos-3, can_ref_pos+1, \
                                  stepper='all', min_mapping_quality=20, \
                                  ignore_overlaps=False, ignore_orphans=False, \
                                  min_base_quality=20):
        if int(mp_line.reference_pos) == can_ref_pos:
            for pileupread in mp_line.pileups:
                if (not pileupread.is_del) and (not pileupread.is_refskip):
                    read_info = pileupread.alignment
                    ref_positions = read_info.get_reference_positions()

                    for q_pos in range(read_info.query_alignment_start, \
                                       read_info.query_alignment_end):
                        var_index = q_pos - read_info.query_alignment_start

                        if var_index < len(ref_positions): # WHY?
                            ref_pos = ref_positions[var_index]
                            ref_sequence = read_info.get_reference_sequence()
                            q_sequence = read_info.query_sequence

                            # Target 위치를 확인하기로 함.
                            if ref_pos == can_ref_pos:
                                if var_index < len(ref_sequence):
                                    ref_base = ref_sequence[var_index]
                                    q_base = q_sequence[q_pos]

                                    # 같은 variant를 가지고 있는 경우임.
                                    if (ref_base.upper() == can_Var_list[0]) and \
                                       (q_base.upper() == can_Var_list[1]):
                                        tmp_var_line = '{0}-{1}'.format(ref_base.upper(), q_base.upper())

                                        # save format: {ref_position: {ref-alt: count}}
                                        if can_ref_pos not in tmp_candidate_dic:
                                            tmp_candidate_dic\
                                                        .setdefault(can_ref_pos, dict())\
                                                        .setdefault(tmp_var_line, 1)
                                        else:
                                            if tmp_var_line not in tmp_candidate_dic[can_ref_pos]:
                                                tmp_candidate_dic[can_ref_pos]\
                                                                 .setdefault(tmp_var_line, 1)
                                            else:
                                                tmp_candidate_dic[can_ref_pos][tmp_var_line] += 1
                                    else: pass
                                else: pass
                            else: pass
                        else: pass
                else: pass
        else: pass

    return tmp_candidate_dic

def non_phase_info_saving(chr_name, start, end, phasing_candidate_Dic, pos_dic, non_phasing_info_Dic):

#    key_line = '{0}_{1}_{2}_{3}'.format(chr_name, pos_info, pos_dic[pos_info]['Ref'], pos_dic[pos_info]['Alt'][0])
    key_line = '{0}_{1}_{2}'.format(chr_name, start, end)

    # Save format: {Chr_name_Pos_info_Ref_Alt: phasing_candidate_Dic}
    non_phasing_info_Dic.setdefault(key_line, phasing_candidate_Dic)

    return non_phasing_info_Dic



def germline_phaser(chr_name, var_list, pos_dic, type_info, c_bam, f_bam, m_bam, o):
    # Load child BAM file.
    samfile_C = pysam.AlignmentFile(c_bam, 'rb', threads=3)

    # Load father BAM file.
    samfile_F = pysam.AlignmentFile(f_bam, 'rb', threads=3)

    # Load mother BAM file.
    samfile_M = pysam.AlignmentFile(m_bam, 'rb', threads=3)

    ####################################################################
    # The phasing results of the DNV will be saved to Python dictionary.
    DNV_phasing_info_Dic = dict(); non_phasing_info_Dic = dict()

    for pos_info in var_list:
        start = int(pos_info.split('-')[0])#
        end = int(pos_info.split('-')[-1])#
        #print(chr_name)
        #print(pos_info)
        #print(pos_dic)
        sv_type = pos_dic[pos_info]['type']#

        extend_flag = 0

        # Check mapped reads in the +-2bp range near the DNV position. # WHY +-2bp range?
        phasing_candidate_Dic = dict(); non_phasing_Dic = dict()

        print('ANALYZING SV START POS')
        for mp_line in samfile_C.pileup(chr_name, start-3, start+1, stepper='all', min_mapping_quality=20, ignore_overlaps=False, ignore_orphans=False, min_base_quality=20):

            # DNV position에 mapping된 read를 식별하기 위한 방법.
            # reference_pos는 0-padded 되었기 때문임.
            if int(mp_line.reference_pos + 1) == start:
                # Collect the mapped reads that have DNM information.
                for pileupread in mp_line.pileups:
                    # For dSNV type
                    #if type_info == 'dSNV':
                    if (not pileupread.is_del) and (not pileupread.is_refskip):
                        read_info = pileupread.alignment
                        ref_positions = read_info.get_reference_positions()
                        ref_sequence = read_info.get_reference_sequence()
                        q_sequence = read_info.query_sequence

                        # 먼저, dSNV를 가진 reads를 찾기로 함.
                        # Find reads being soft-clipped on the SV start / end position
                        clip_tag = 0
                        if is_clipped(read_info, start, end) or is_discordant(read_info, chr_name, start, end, sv_type, samfile_C):#
                            clip_tag += 1;

#                        for q_pos in range(read_info.query_alignment_start, read_info.query_alignment_end):
#                            var_index = q_pos - read_info.query_alignment_start
#
#                            if var_index < len(ref_positions):
#                                ref_pos = ref_positions[var_index]
#
#                                # DNV 정보를 가짐을 확인함.
#                                if int(ref_pos+1) == start:
#                                    if var_index < len(ref_sequence):
#                                        ref_base = ref_sequence[var_index]
#                                        q_base = q_sequence[q_pos]
#
#                                        if q_base.upper() == pos_dic[pos_info]['Alt'][0].upper():
#                                            DNV_tag += 1; break
#                                        else: pass
#                                    else: pass
#
#                                # Pass the non-target reads
#                                else: pass
#                            else: pass

                        # DNV를 가진 read만 확인함.
                        # 이때, mate-read에서도 추가적인 정보를 찾아봐야 할 것 같음.
                        if clip_tag > 0:
                            # 해당 read 내에 informative site 정보를 모두 수집하기로 함.
                            phasing_candidate_Dic = candidate_Vars_parser(read_info, start, phasing_candidate_Dic)

                            # 해당 read의 pair 정보를 가져옴.
                            if read_info.is_paired:

                                # 같은 chromosome 내에 mapping이 되어야 같은 phase 정보를 얻을 수 있지 않을까?
                                if read_info.next_reference_name == chr_name:

                                    # 해당 read의 mate-read 정보를 확인하기로 함.
                                    phasing_candidate_Dic = fetch_mate(start, c_bam, \
                                                                       read_info, phasing_candidate_Dic)

                                else: pass
                            else: pass

                        # DNV를 가진 read가 아닌 경우임.
                        else:
                            # 해당 read 내에 informative site 정보를 모두 수집하기로 함.
                            non_phasing_Dic = candidate_Vars_parser(read_info, start, non_phasing_Dic)

                            # 해당 read의 pair 정보를 가져옴.
                            if read_info.is_paired:

                                # 같은 chromosome 내에 mapping이 되어야 같은 phase 정보를 얻을 수 있지 않을까?
                                if read_info.next_reference_name == chr_name:

                                    # 해당 read의 mate-read 정보를 확인하기로 함.
                                    non_phasing_Dic = fetch_mate(start, c_bam, \
                                                                 read_info, non_phasing_Dic)

                                else: pass
                            else: pass
                        #else: pass

#                    # For dINDEL type
#                    else:
#                        # 먼저, dINDEL를 가진 reads를 찾기로 함.
#                        if pileupread.indel == len(pos_dic[pos_info]['Alt'][0])-len(pos_dic[pos_info]['Ref']):
#                            read_info = pileupread.alignment
#
#                            # 해당 read 내에 informative site 정보를 모두 수집하기로 함.
#                            phasing_candidate_Dic = candidate_Vars_parser(read_info, pos_info, \
#                                                                          phasing_candidate_Dic)
#
#                            # 해당 read의 pair 정보를 가져옴.
#                            if read_info.is_paired:
#
#                                # 같은 chromosome 내에 mapping이 되어야 같은 phase 정보를 얻을 수 있지 않을까?
#                                if read_info.next_reference_name == chr_name:
#
#                                    # 해당 read의 mate-read 정보를 확인하기로 함.
#                                    phasing_candidate_Dic = fetch_mate(pos_info, c_bam, \
#                                                                       read_info, phasing_candidate_Dic)
#                                else: pass
#                            else: pass
#
#                        # DNV를 가지지 않은 read가 아님.
#                        else:
#                            read_info = pileupread.alignment
#
#                            # 해당 read 내에 informative site 정보를 모두 수집하기로 함.
#                            non_phasing_Dic = candidate_Vars_parser(read_info, pos_info, non_phasing_Dic)
#
#                            # 해당 read의 pair 정보를 가져옴.
#                            if read_info.is_paired:
#
#                                # 같은 chromosome 내에 mapping이 되어야 같은 phase 정보를 얻을 수 있지 않을까?
#                                if read_info.next_reference_name == chr_name:
#
#                                    # 해당 read의 mate-read 정보를 확인하기로 함.
#                                    non_phasing_Dic = fetch_mate(pos_info, c_bam, \
#                                                                 read_info, non_phasing_Dic)
#
#                                else: pass
#                            else: pass
            else: pass
        print('ANALYZING SV END POS')
        for mp_line in samfile_C.pileup(chr_name, end-3, end+1, \
                                        stepper='all', min_mapping_quality=20, \
                                        ignore_overlaps=False, ignore_orphans=False, \
                                        min_base_quality=20):

            # DNV position에 mapping된 read를 식별하기 위한 방법.
            # reference_pos는 0-padded 되었기 때문임.
            if int(mp_line.reference_pos + 1) == end:
                # Collect the mapped reads that have DNM information.
                for pileupread in mp_line.pileups:
                    # For dSNV type
                    #if type_info == 'dSNV':
                    if (not pileupread.is_del) and (not pileupread.is_refskip):
                        read_info = pileupread.alignment
                        ref_positions = read_info.get_reference_positions()
                        ref_sequence = read_info.get_reference_sequence()
                        q_sequence = read_info.query_sequence

                        # 먼저, dSNV를 가진 reads를 찾기로 함.
                        # Find reads being soft-clipped on the SV start / end position
                        clip_tag = 0
                        if is_clipped(read_info, end, start) or is_discordant(read_info, chr_name, end, start, sv_type, samfile_C):#
                            clip_tag += 1

                        # DNV를 가진 read만 확인함.
                        # 이때, mate-read에서도 추가적인 정보를 찾아봐야 할 것 같음.
                        if clip_tag > 0:
                            # 해당 read 내에 informative site 정보를 모두 수집하기로 함.
                            phasing_candidate_Dic = candidate_Vars_parser(read_info, end, phasing_candidate_Dic)

                            # 해당 read의 pair 정보를 가져옴.
                            if read_info.is_paired:

                                # 같은 chromosome 내에 mapping이 되어야 같은 phase 정보를 얻을 수 있지 않을까?
                                if read_info.next_reference_name == chr_name:

                                    # 해당 read의 mate-read 정보를 확인하기로 함.
                                    phasing_candidate_Dic = fetch_mate(end, c_bam, \
                                                                       read_info, phasing_candidate_Dic)

                                else: pass
                            else: pass

                        # DNV를 가진 read가 아닌 경우임.
                        else:
                            # 해당 read 내에 informative site 정보를 모두 수집하기로 함.
                            non_phasing_Dic = candidate_Vars_parser(read_info, end, non_phasing_Dic)

                            # 해당 read의 pair 정보를 가져옴.
                            if read_info.is_paired:

                                # 같은 chromosome 내에 mapping이 되어야 같은 phase 정보를 얻을 수 있지 않을까?
                                if read_info.next_reference_name == chr_name:

                                    # 해당 read의 mate-read 정보를 확인하기로 함.
                                    non_phasing_Dic = fetch_mate(end, c_bam, \
                                                                 read_info, non_phasing_Dic)

                                else: pass
                            else: pass


        # 여기서, DNV read만 가지고 있는 variant 정보만 남김.
        phasing_candidate_Dic = check_only_DNV_reads(phasing_candidate_Dic, non_phasing_Dic)

        # DNV 가진 read 내에 phasing을 위한 후보군 변이 정보가 한 개라도 있었던 경우.
        # 해당 정보에서 부모와의 phasing 가능한 경우가 있는지 확인하기로 함.
        phasing_tag_Dic = dict()

        if len(phasing_candidate_Dic) > 0:
            # 결국 하나의 DNV에 여러 개의 position이 나올 것이기 때문에,
            # 모든 정보를 모은 다음에, 일관성 있게 정보가 나왔는지 확인하기로 함.
            # Save format: {ref_position: {Ref-Alt: Count}}
            for can_ref_pos, can_Var_Dic in phasing_candidate_Dic.items():
                # DNV를 가진 read들이 해당 위치에 단일 변이를 가지고 있음.
                # 해당 위치에 여러 종류의 변이가 존재하면, error일 확률이 높음.
                if len(can_Var_Dic) == 1:
                    can_Var_list = [a for a in can_Var_Dic.keys()][0].split('-') # Ref-Alt

                    # 아버지 쪽 확인하기
                    f_candidate_dic = collect_the_variants_in_specific_position(samfile_F, chr_name, can_ref_pos, can_Var_list)
                    # 어머니쪽 확인하기
                    m_candidate_dic = collect_the_variants_in_specific_position(samfile_M, chr_name, can_ref_pos, can_Var_list)

                    # 해당 위치의 변이가 부모 중 한 쪽에서만 발견된다면?
                    # 그건 DNV가 해당 부모에서 유래한 것으로 유추할 수 있지 않을까?
                    # Save format: {ref_position: {ref-alt: count}}
                    # Save format: {F: {ref_position: {ref-alt: count}}}
                    if len(f_candidate_dic) > 0:
                        if 'Father' not in phasing_tag_Dic:
                            phasing_tag_Dic.setdefault('Father', f_candidate_dic)
                        else:
                            phasing_tag_Dic['Father'].update(f_candidate_dic)
                    else: pass

                    # Save format: {M: {ref_position: {ref-alt: count}}}
                    if len(m_candidate_dic) > 0:
                        if 'Mother' not in phasing_tag_Dic:
                            phasing_tag_Dic.setdefault('Mother', m_candidate_dic)
                        else:
                            phasing_tag_Dic['Mother'].update(m_candidate_dic)
                    else: pass

                # 해당 위치에 여러 개의 변이 정보가 있으면 확인이 안되는 케이스로 분류.
                else: pass

#################### checked 240828

        # Extended mode를 사용해서 다시 한 번 찾아보기로 함.
        else: extend_flag += 1

        # 해당 DNV를 가진 read에서 발견된 후보군 위치 정보들을 모두 확인했을 때,
        # 아버지 또는 어머니 중 한 명에서만 정보가 link되었기 때문에
        # 해당 DNV가 link된 부모 쪽에서 유래되었을 가능성이 생김.
        if len(phasing_tag_Dic) == 1:
            # Save format: {F or M: {ref_position: {ref-alt: count}}}
            for phasing_info, phasing_var_Dic in phasing_tag_Dic.items():
                # Save format: {F or M: {Chr_name: {Pos_info: {ref_position: {ref-alt: count}}}}}
                if phasing_info not in DNV_phasing_info_Dic:
                    DNV_phasing_info_Dic.setdefault(phasing_info, dict()).setdefault(chr_name, dict()).setdefault(start, phasing_var_Dic)
                else:
                    if chr_name not in DNV_phasing_info_Dic[phasing_info]:
                        DNV_phasing_info_Dic[phasing_info].setdefault(chr_name, dict()).setdefault(start, phasing_var_Dic)
                    else:
                        DNV_phasing_info_Dic[phasing_info][chr_name].setdefault(start, phasing_var_Dic)

        # Father or Mother에서 확인해본 결과 DNV read 내 variant들이 아버지와 어머니
        # 모두에게서 하나 이상 발견되었음.
        elif len(phasing_tag_Dic) > 1:
            phasing_stage_info = 'None'

            # 아버지에서 해당 위치의 variant들을 1개 초과로 발견하면서, 어머니 쪽에서는
            # 한 개 이하로 발견된 상황이면, Father 쪽으로 phasing함.
            #if len(phasing_tag_Dic['Father']) > 1 and len(phasing_tag_Dic['Mother']) <= 1:
            #if len(phasing_tag_Dic['Father']) > len(phasing_tag_Dic['Mother']):
            if (len(phasing_tag_Dic['Father']) > len(phasing_tag_Dic['Mother'])) and (len(phasing_tag_Dic['Mother']) <= 1):
                phasing_stage_info = 'Father'

            # 어머니에서 해당 위치의 variant들을 1개 초과로 발견하면서, 아버지 쪽에서는
            # 한 개 이하로 발견된 상황이면, Mother 쪽으로 phasing함.
            #elif len(phasing_tag_Dic['Father']) < len(phasing_tag_Dic['Mother']):
            if (len(phasing_tag_Dic['Mother']) > len(phasing_tag_Dic['Father'])) and (len(phasing_tag_Dic['Father']) <= 1):
                phasing_stage_info = 'Mother'

            # 이도 저도 아닌 상황임.
            else: pass
                
            if phasing_stage_info != 'None':
                # Save format: {F or M: {Chr_name: {Pos_info: {ref_position: {ref-alt: count}}}}}
                if phasing_stage_info not in DNV_phasing_info_Dic:
                    DNV_phasing_info_Dic.setdefault(phasing_stage_info, dict()).setdefault(chr_name, dict()).setdefault(start, phasing_tag_Dic[phasing_stage_info])
                else:
                    if chr_name not in DNV_phasing_info_Dic[phasing_stage_info]:
                        DNV_phasing_info_Dic[phasing_stage_info].setdefault(chr_name, dict()).setdefault(start, phasing_tag_Dic[phasing_stage_info])
                    else:
                        DNV_phasing_info_Dic[phasing_stage_info][chr_name].setdefault(start, phasing_tag_Dic[phasing_stage_info])

            else: extend_flag += 1
        else: extend_flag += 1

        # Extended mode 사용을 위해서, DNV read에서만 발견되는 variant 정보를 저장하기로 함.
        if extend_flag > 0:
            # Save format: {Chr_name_Pos_info_Ref_Alt: phasing_candidate_Dic}
            non_phasing_info_Dic = non_phase_info_saving(chr_name, start, end, phasing_candidate_Dic, pos_dic, non_phasing_info_Dic)

        # Read extension mode를 사용하지 않아도 되는 경우임.
        else: pass

    o.put([DNV_phasing_info_Dic, non_phasing_info_Dic])

    samfile_C.close()
    samfile_F.close()
    samfile_M.close()


    # NEED TO ADD CODES FOR SV END POSITION
    ## ADDED ON LINE 361 


#################### checked 240829
def extension_info_saver(c_bam, extension_info_Dic, DNV_info_list, phasing_candidate_Dic):
    # Load the sample bam files
    samfile_C = pysam.AlignmentFile(c_bam, 'rb', threads=2)

    # 해당 DNV read만 가지는 variant 정보임.
    for tmp_pos, tmp_var_Dic in extension_info_Dic.items(): # extension_info_Dic in the form of {18032628: {'A-T': 1}, 18032650: {'T-A': 1}}
        tmp_var_list = [z.split('-') for z in tmp_var_Dic.keys()][0]

        for mp_line in samfile_C.pileup(DNV_info_list[0], int(tmp_pos)-3, int(tmp_pos)+1, stepper='all', min_mapping_quality=20, ignore_overlaps=False, ignore_orphans=False, min_base_quality=20):

            if mp_line.reference_pos == int(tmp_pos):
                read_limit_flag = 0
                for pileupread in mp_line.pileups:
                    if read_limit_flag <= 1000:
                        read_limit_flag += 1

                        if (not pileupread.is_del) and (not pileupread.is_refskip):
                            read_info = pileupread.alignment
                            ref_positions = read_info.get_reference_positions()
                            ref_sequence = read_info.get_reference_sequence()
                            q_sequence = read_info.query_sequence

                            for q_pos in range(read_info.query_alignment_start, read_info.query_alignment_end):
                                var_index = q_pos - read_info.query_alignment_start

                                if var_index < len(ref_positions):
                                    ref_pos = ref_positions[var_index]

                                    if ref_pos == tmp_pos:
                                        if var_index < len(ref_sequence):
                                            ref_base = ref_sequence[var_index]
                                            q_base = q_sequence[q_pos]

                                            # DNV read가 가지고 있는 variant 가지고 있음.
                                            if (ref_base.upper() == tmp_var_list[0]) and \
                                               (q_base.upper() == tmp_var_list[1]):

                                                # 해당 read 내에 informative site 정보를 모두 수집하기로 함.
                                                phasing_candidate_Dic = candidate_Vars_parser(read_info, int(DNV_info_list[1]), \
                                                                                              phasing_candidate_Dic)

                                                # 해당 read의 pair 정보를 가져옴.
                                                if read_info.is_paired:

                                                    # 같은 chromosome 내에 mapping이 되어야 같은 phase 정보를 얻을 수 있지 않을까?
                                                    if read_info.next_reference_name == DNV_info_list[0]:

                                                        # 해당 read의 mate-read 정보를 확인하기로 함.
                                                        phasing_candidate_Dic = fetch_mate(int(DNV_info_list[1]), c_bam, \
                                                                                           read_info, phasing_candidate_Dic)

                                                    else: pass
                                                else: pass
                                            else: pass
                                        else: pass
                                    else: pass
                                else: pass
                        else: pass
                    else: break
            else: pass

    samfile_C.close()

    return phasing_candidate_Dic


def read_extension_phaser(check_vDic, for_Vars_list, non_phasing_info_Dic, \
                             c_bam, f_bam, m_bam, cpu_num, o):

    # 아버지와 어머니 bam file을 load하자!
    samfile_F = pysam.AlignmentFile(f_bam, 'rb', threads=2)

    samfile_M = pysam.AlignmentFile(m_bam, 'rb', threads=2)

    DNV_phasing_info_Dic = dict()

    for non_var_line in for_Vars_list:
        # DNV 정보임.
        # Format: Chr_name_PosInfo_Ref_Alt
        non_var_list = non_var_line.split('_')

        none_flag = 0; phasing_candidate_Dic = dict()

        # 해당 DNV read만 가지는 variant 정보임.
        # 만약 한 개라도 그 정보를 가지고 있지 않다면, 불가능하겠지?
        if len(non_phasing_info_Dic[non_var_line]) > 0:
            phasing_candidate_Dic = extension_info_saver(c_bam, non_phasing_info_Dic[non_var_line], \
                                                         non_var_list, phasing_candidate_Dic)

        else: pass # DNV read가 자신만의 variant를 가지지 못했음.

        # 한 번 더 extension 시켜보자.
        if len(phasing_candidate_Dic) > 0:
            new_phasing_candidate_Dic = dict()

            # 한 번 더 extension 시켜서 정보를 가져오기로 함.
            new_phasing_candidate_Dic = extension_info_saver(c_bam, phasing_candidate_Dic, \
                                                             non_var_list, new_phasing_candidate_Dic)

            if len(new_phasing_candidate_Dic) > 0:
                for tmp_new_pos, tmp_new_var_Dic in new_phasing_candidate_Dic.items():
                    if tmp_new_pos not in phasing_candidate_Dic:
                        phasing_candidate_Dic.setdefault(tmp_new_pos, tmp_new_var_Dic)
                    else:
                        for tmp_new_var_line, tmp_new_var_count in tmp_new_var_Dic.items():
                            if tmp_new_var_line not in phasing_candidate_Dic[tmp_new_pos]:
                                phasing_candidate_Dic[tmp_new_pos].setdefault(tmp_new_var_line, tmp_new_var_count)
                            else:
                                phasing_candidate_Dic[tmp_new_pos][tmp_new_var_line] += tmp_new_var_count
            else: pass

        # extension 시켰을 때, variant가 발견되지 않음.
        else: pass

        # DNV 가진 read 내에 phasing을 위한 후보군 변이 정보가 한 개라도 있었던 경우.
        # 해당 정보에서 부모와의 phasing 가능한 경우가 있는지 확인하기로 함.
        phasing_tag_Dic = dict()

        if len(phasing_candidate_Dic) > 0:
            # Save format: {ref_position: {Ref-Alt: Count}}
            for can_ref_pos, can_Var_Dic in phasing_candidate_Dic.items():
                if len(can_Var_Dic) == 1:
                    can_Var_list = [a for a in can_Var_Dic.keys()][0].split('-') # Ref-Alt

                    # 아버지 쪽 확인하기
                    f_candidate_dic = collect_the_variants_in_specific_position(\
                                            samfile_F, non_var_list[0], can_ref_pos, can_Var_list)
                    # 어머니쪽 확인하기
                    m_candidate_dic = collect_the_variants_in_specific_position(\
                                            samfile_M, non_var_list[0], can_ref_pos, can_Var_list)

                    # Save format: {ref_position: {ref-alt: count}}
                    # Save format: {F: {ref_position: {ref-alt: count}}}
                    if len(f_candidate_dic) > 0:
                        # Save format: {F: {ref_position: {ref-alt: count}}}
                        if 'Father' not in phasing_tag_Dic:
                            phasing_tag_Dic.setdefault('Father', f_candidate_dic)
                        else:
                            phasing_tag_Dic['Father'].update(f_candidate_dic)
                    else: pass

                    if len(m_candidate_dic) > 0:
                        # Save format: {M: {ref_position: {ref-alt: count}}}
                        if 'Mother' not in phasing_tag_Dic:
                            phasing_tag_Dic.setdefault('Mother', m_candidate_dic)
                        else:
                            phasing_tag_Dic['Mother'].update(m_candidate_dic)
                    else: pass

                else: pass

        # Read extension mode를 사용해도 찾을 수 없는 상황임.
        else: none_flag += 1

        if len(phasing_tag_Dic) == 1:
            # Save format: {F or M: {ref_position: {ref-alt: count}}}
            for phasing_info, phasing_var_Dic in phasing_tag_Dic.items():
                # Save format: {F or M: {Chr_name: {Pos_info: {ref_position: {ref-alt: count}}}}}
                if phasing_info not in DNV_phasing_info_Dic:
                    DNV_phasing_info_Dic.setdefault(phasing_info, dict())\
                                        .setdefault(non_var_list[0], dict())\
                                        .setdefault(int(non_var_list[1]), phasing_var_Dic)
                else:
                    if non_var_list[0] not in DNV_phasing_info_Dic[phasing_info]:
                        DNV_phasing_info_Dic[phasing_info]\
                                        .setdefault(non_var_list[0], dict())\
                                        .setdefault(int(non_var_list[1]), phasing_var_Dic)
                    else:
                        DNV_phasing_info_Dic[phasing_info][non_var_list[0]]\
                                        .setdefault(int(non_var_list[1]), phasing_var_Dic)

        elif len(phasing_tag_Dic) > 1:
            phasing_stage_info = 'None'

            if len(phasing_tag_Dic['Father']) > len(phasing_tag_Dic['Mother']):
                phasing_stage_info = 'Father'

            elif len(phasing_tag_Dic['Father']) < len(phasing_tag_Dic['Mother']):
                phasing_stage_info = 'Mother'

            else: pass
                
            if phasing_stage_info != 'None':
                # Save format: {F or M: {Chr_name: {Pos_info: {ref_position: {ref-alt: count}}}}}
                if phasing_stage_info not in DNV_phasing_info_Dic:
                    DNV_phasing_info_Dic.setdefault(phasing_stage_info, dict())\
                                        .setdefault(non_var_list[0], dict())\
                                        .setdefault(int(non_var_list[1]), phasing_tag_Dic[phasing_stage_info])
                else:
                    if non_var_list[0] not in DNV_phasing_info_Dic[phasing_stage_info]:
                        DNV_phasing_info_Dic[phasing_stage_info]\
                                        .setdefault(non_var_list[0], dict())\
                                        .setdefault(int(non_var_list[1]), phasing_tag_Dic[phasing_stage_info])
                    else:
                        DNV_phasing_info_Dic[phasing_stage_info][non_var_list[0]]\
                                        .setdefault(int(non_var_list[1]), phasing_tag_Dic[phasing_stage_info])

            else: none_flag += 1
        else: none_flag += 1

        # Extension mode를 사용해도 여전히 parent-of-origin을 찾을 수 없음.
        if none_flag > 0:
            # Save format: {None: {Chr_name: {Pos_info: {ref_position: {ref-alt: count}}}}}
            if 'None' not in DNV_phasing_info_Dic:
                DNV_phasing_info_Dic.setdefault('None', dict())\
                                    .setdefault(non_var_list[0], dict())\
                                    .setdefault(int(non_var_list[1]), dict())
            else:
                if non_var_list[0] not in DNV_phasing_info_Dic['None']:
                    DNV_phasing_info_Dic['None']\
                                    .setdefault(non_var_list[0], dict())\
                                    .setdefault(int(non_var_list[1]), dict())
                else:
                    DNV_phasing_info_Dic['None'][non_var_list[0]]\
                                    .setdefault(int(non_var_list[1]), dict())
        else: pass

    o.put(DNV_phasing_info_Dic)

def phase_by_germline_info(check_vDic, type_info, c_bam, f_bam, m_bam, \
                             cpu_num, var_count):
#def phase_by_germline_info(v_l, type_info, c_bam, f_bam, m_bam, cpu_num, var_count): 

    #phasing results of the DNV will be saved to Python dictionary.
    DNV_phasing_info_Dic = dict(); non_phasing_info_Dic = dict()

    for chr_name, pos_dic in check_vDic.items():# check_vDic format = {contig : {start_pos-end_pos : {sv_type : SVTYPE}}}
        p_list = list(); o_list = list()

        # The Case of variant positions less than cpu number
        if len(pos_dic) <= cpu_num:

            for e_pos in pos_dic.keys():
                o = Queue(); o_list.append(o)

                p_list.append(Process(target = germline_phaser, \
                              args = (chr_name, [e_pos], pos_dic, type_info, \
                                      c_bam, f_bam, m_bam, o,)))

        # The Case of variant positions more than cpu number
        else:
            split_num = math.floor(len(pos_dic)/float(cpu_num))
            var_list = [x for x in pos_dic.keys()]; step_num = 0

            for a in range(0, cpu_num+1):
                o = Queue(); o_list.append(o)

                if a == 0:
                    p_list.append(Process(target = germline_phaser, \
                                  args = (chr_name, var_list[0:split_num], pos_dic, type_info, \
                                          c_bam, f_bam, m_bam, o,)))
                    step_num += split_num
                elif a < cpu_num:
                    p_list.append(Process(target = germline_phaser, \
                                  args = (chr_name, var_list[step_num:step_num+split_num], pos_dic, type_info, \
                                          c_bam, f_bam, m_bam, o,)))
                    step_num += split_num
                else:
                    p_list.append(Process(target = germline_phaser, \
                                  args = (chr_name, var_list[step_num:], pos_dic, type_info, \
                                          c_bam, f_bam, m_bam, o,)))

        for each_p in p_list: each_p.start()

        for each_o in o_list:
            r = each_o.get()

            # Save format: {Father or Mother: {Chr_name: {Pos_info: {ref_position: {ref-alt: count}}}}}
            if len(r[0]) > 0:
                for phasing_info, phasing_Var_dic in r[0].items():
                    if phasing_info not in DNV_phasing_info_Dic:
                        DNV_phasing_info_Dic\
                                        .setdefault(phasing_info, phasing_Var_dic)
                    else:
                        for phasing_chr, phasing_Var_info_dic in phasing_Var_dic.items():
                            if phasing_chr not in DNV_phasing_info_Dic[phasing_info]:
                                DNV_phasing_info_Dic[phasing_info]\
                                                .setdefault(phasing_chr, phasing_Var_info_dic)
                            else:
                                for phasing_pos, phasing_detail_info_dic in phasing_Var_info_dic.items():
                                    DNV_phasing_info_Dic[phasing_info][phasing_chr]\
                                                    .setdefault(phasing_pos, phasing_detail_info_dic)
            else: pass

            # Save format: {Chr_name_Pos_info_Ref_Alt: phasing_candidate_Dic}
            if len(r[1]) > 0:
                for re_line, re_Var_Dic in r[1].items():
                    # Save format: {Chr_name_Pos_info: {'Ref':A}}}
                    non_phasing_info_Dic\
                                .setdefault(re_line, re_Var_Dic)
            else: pass

            del r; each_o.close()

        for each_p in p_list: each_p.join()

        del p_list; del o_list

    #print ('# This {0} DNV has parental phase information.'\
    #                    .format(var_count - len(non_phasing_info_Dic)))

    if 'Father' in DNV_phasing_info_Dic:
        tmp_c = 0
        for x, y in DNV_phasing_info_Dic['Father'].items():
            print (x, len(y))
            tmp_c += len(y)
        print ('Father:{0}'.format(tmp_c))
    if 'Mother' in DNV_phasing_info_Dic:
        tmp_c = 0
        for x, y in DNV_phasing_info_Dic['Mother'].items():
            print (x, len(y))
            tmp_c += len(y)
        print ('Mother:{0}'.format(tmp_c))
    if 'None' in DNV_phasing_info_Dic:
        tmp_c = 0
        for x, y in DNV_phasing_info_Dic['None'].items():
            print (x, len(y))
            tmp_c += len(y)
        print ('None:{0}'.format(tmp_c))

    ###########################################
    # DNV 중 phasing되지 않은 것들이 있음.
    # 이들은 read extension mode를 적용해야 함.
    if len(non_phasing_info_Dic) > 0:
        print(non_phasing_info_Dic)
        print ('#  {0} DNVs must use the read extension mode of the DNV phase estimator.'\
                                                           .format(len(non_phasing_info_Dic)))

        # 이번 전략은 DNV에서만 발견되었던 variant 정보를 가진 read를 가져오는 방식으로
        # 이어 붙이기를 해보고, 그 read에서 variant 정보를 추가로 가져오는 방식을 선택하기로 함.
        p_list = list(); o_list = list(); new_non_phasing_info_Dic = dict()

        if len(non_phasing_info_Dic) > cpu_num:
            split_flag = math.floor(len(non_phasing_info_Dic)/float(cpu_num))
            step_index = 0; for_non_Vars_list = [a for a in non_phasing_info_Dic.keys()]

            for x in range(0, cpu_num+1):
                o = Queue(); o_list.append(o)

                if x == 0:
                    p_list.append(Process(target = read_extension_phaser, \
                                          args = (check_vDic, for_non_Vars_list[0:split_flag], \
                                                  non_phasing_info_Dic, \
                                                  c_bam, f_bam, m_bam, cpu_num, o,)))
                    step_index += split_flag
                elif x < cpu_num:
                    p_list.append(Process(target = read_extension_phaser, \
                                          args = (check_vDic, for_non_Vars_list[step_index:step_index+split_flag], \
                                                  non_phasing_info_Dic, \
                                                  c_bam, f_bam, m_bam, cpu_num, o,)))
                    step_index += split_flag
                else:
                    p_list.append(Process(target = read_extension_phaser, \
                                          args = (check_vDic, for_non_Vars_list[step_index:], \
                                                  non_phasing_info_Dic, \
                                                  c_bam, f_bam, m_bam, cpu_num, o,)))
        else:
            for non_var_line in non_phasing_info_Dic.keys():
                o = Queue(); o_list.append(o)
                p_list.append(Process(target = read_extension_phaser, \
                                      args = (check_vDic, [non_var_line], \
                                              non_phasing_info_Dic, \
                                              c_bam, f_bam, m_bam, cpu_num, o,)))

        for e_p in p_list: e_p.start()

        for each_o in o_list:
            r = each_o.get()

            # Save format: {Father or Mother: {Chr_name: {Pos_info: {ref_position: {ref-alt: count}}}}}
            if len(r) > 0:
                for phasing_info, phasing_Var_dic in r.items():
                    if phasing_info not in DNV_phasing_info_Dic:
                        DNV_phasing_info_Dic\
                                        .setdefault(phasing_info, phasing_Var_dic)
                    else:
                        for phasing_chr, phasing_Var_info_dic in phasing_Var_dic.items():
                            if phasing_chr not in DNV_phasing_info_Dic[phasing_info]:
                                DNV_phasing_info_Dic[phasing_info]\
                                                .setdefault(phasing_chr, phasing_Var_info_dic)
                            else:
                                for phasing_pos, phasing_detail_info_dic in phasing_Var_info_dic.items():
                                    DNV_phasing_info_Dic[phasing_info][phasing_chr]\
                                                    .setdefault(phasing_pos, phasing_detail_info_dic)
            else: pass

            del r; each_o.close()

        for each_p in p_list: each_p.join()

    # 모두 phasing을 했다고???
    else: pass

    if 'Father' in DNV_phasing_info_Dic:
        tmp_c = 0
        for x, y in DNV_phasing_info_Dic['Father'].items():
            print (x, len(y))
            tmp_c += len(y)
        print ('Father:{0}'.format(tmp_c))
    if 'Mother' in DNV_phasing_info_Dic:
        tmp_c = 0
        for x, y in DNV_phasing_info_Dic['Mother'].items():
            print (x, len(y))
            tmp_c += len(y)
        print ('Mother:{0}'.format(tmp_c))
    if 'None' in DNV_phasing_info_Dic:
        tmp_c = 0
        for x, y in DNV_phasing_info_Dic['None'].items():
            print (x, len(y))
            tmp_c += len(y)
        print ('None:{0}'.format(tmp_c))

        print ('# {0} DNVs were phased by the read extension mode of the DNV phase estimator.'\
                                                           .format(len(non_phasing_info_Dic)-tmp_c))

    return DNV_phasing_info_Dic

def result_printer(phased_d, type_info, outname):    
    outfile = open('{0}.parent_of_origin.table'.format(outname), 'w')
    outfile.write('#Phasing\tSV_type\tcontig\tstart\tPhasing_Var_position\tPhasing_Var_info\tPhasing_Var_count\n')

    # Save format: {F or M: {contig: {Pos_info: {ref_position: {ref-alt: count}}}}}
    for phasing_info, phasing_pos_d in phased_d.items():
        for DNV_chr, phasing_var_d in phasing_pos_d.items():
            for DNV_pos, phasing_vars_d in phasing_var_d.items():
                tmp_out_list = list()
                for phasing_pos, phasing_detail_d in phasing_vars_d.items():
                    for phasing_var, phasing_count in phasing_detail_d.items():
                        tmp_out_list.append(str(phasing_pos+1))
                        tmp_out_list.append(phasing_var)
                        tmp_out_list.append(str(phasing_count))
                #print(DNV_pos)
                outfile.write('{0}\t{1}\t{2}\t{3}\n'.format(\
                              phasing_info, DNV_chr, DNV_pos, '\t'.join(tmp_out_list)))

    outfile.close()

def is_clipped(read, start, end):
    """
    | A FUNCTION THAT RETURNS BOOLEAN WHETHER A GIVEN READ IS EXACTLY CLIPPED ON THE SV START/END POSITION
    : read : pysam read object
    : start : SV start position
    : end : SV end position
    | RETURNS : True / False |
    """

    result = False

    if ((not read.is_unmapped) and (read.mapping_quality >= 20) and (read.mapping_quality != 255) \
        and (not read.flag & 0x100) \
        and (not read.flag & 0x400)): # 0x100 : not primary-alignment, 0x400 : read is PCR or opticalduplicate

            cigar = read.cigarstring
            align_start = read.reference_start
            align_end = read.reference_end -1
            if 'S' in cigar:
                #print(start, cigar)
                cigar_match = re.findall(r'(\d+)([A-Z])', cigar)
                clip_findex = cigar.find('S')
                clip_rindex = cigar.rfind('S')
                match_index = cigar.find('M')

                n_match = [int(match[0]) for match in cigar_match if match[1] not in ('S', 'I')]
                len_match = sum(n_match)

                if (clip_findex > match_index) and (clip_rindex > match_index):
                    clip_pos = align_start + len_match -1

                elif (clip_findex < match_index) and (clip_rindex > match_index):
                    if start < end:
                        clip_pos = align_start
                    else:
                        clip_pos = align_start + len_match -1

                elif (clip_findex < match_index) and (clip_rindex < match_index):
                    clip_pos = align_start
                
                #print('clip_read_function', align_start, clip_pos, start)
                #if clip_pos == start:
                if (clip_pos >= (start - 1)) and (clip_pos <= (start + 1)):
                    result = True

    return result

def is_discordant(read, chr_name, start, end, sv_type, bamobj):
    """
    | A FUNCTION THAT RETURNS BOOLEAN WHETHER A GIVEN READ IS A DISCORDANT READ PAIR SPANNING THE SV
    : read : pysam read object
    : start : SV start position
    : end : SV end position
    | RETURNS : True / False 
    """

    result = False
    discordant = False; mate_match = False
    breakend_extend_size = 300
#    if read.next_reference_name == mate_contig:
    if read.next_reference_name == chr_name:
        if (read.flag & 0x1) and (not read.flag & 0x4) \
        and (not read.flag & 0x8) and (not read.flag & 0x100) \
        and (not read.flag & 0x400) and ((not read.flag & 0x2) \
        or ((read.flag & 0x10 and read.flag & 0x20) or ((not read.flag & 0x10) \
        and (not read.flag & 0x20))) or (read.next_reference_start < read.reference_start)):
        #if ((read.is_paired) and (not read.mate_is_unmapped) and (not read.is_proper_pair)) or \
        #((read.is_paired) and (not read.mate_is_unmapped) and (read.is_reverse == read.mate_is_reverse)):
            if sv_type == 'INV':
                #if read.is_reverse == read.mate_is_reverse:
                if (read.flag & 0x10 and read.flag & 0x20) or \
                ((not read.flag & 0x10) and (not read.flag & 0x20)):
                    discordant_bool = True
                    discordant_type = 'orient'

            elif sv_type == 'DEL':
                if (not read.flag & 0x2) or \
                (abs(read.template_length) > 850):
                    discordant_bool = True
                    discordant_type = 'size'

            elif sv_type == 'DUP':
                rl_pair = False

                try:
                    mate_read = bamobj.mate(read)
                except Exception as e:
                    print('an exception occured, ', e)
                    if read.flag & 0x40: # first in pair
                        if (read.next_reference_start < read.reference_start):
                            rl_pair = True
                    elif read.flag & 0x80: # second in pair
                        if (read.reference_start < read.next_reference_start):
                            rl_pair = True
                try:
                    mate_read = bamobj.mate(read)
                except Exception as e:
                    print('an exception occured, ', e)
                    if read.flag & 0x40: # first in pair
                        if (read.next_reference_start < read.reference_start):
                            rl_pair = True
                    elif read.flag & 0x80: # second in pair
                        if (read.reference_start < read.next_reference_start):
                            rl_pair = True
                else:
                    if read.flag & 0x40: # first in pair
                        if (mate_read.reference_start < read.reference_start):
                            rl_pair = True
                    elif read.flag & 0x80: # second in pair
                        if (read.reference_start < mate_read.reference_start):
                            rl_pair = True
                    else:
                        print(read.query_name, read.flag)
                finally:
                    if rl_pair or \
                    (not read.flag & 0x2) or \
                    (abs(read.template_length) > 850) or (abs(read.template_length) < 250):
                        discordant_bool = True
                        discordant_type = 'size'
            try:
                mate_read = bamobj.mate(read)
            except Exception as e:
                if (read.next_reference_start >= end - breakend_extend_size) and \
                (read.next_reference_end <= end + breakend_extend_size):
                    bp_match = True

            else:
                if (mate_read.reference_start >= end - breakend_extend_size) and \
                (mate_read.reference_end <= end + breakend_extend_size):
                    bp_match = True
            finally:
                if discordant and bp_match:
                    result = True

    return result
 

def main(args):
    check_vDic, var_count = variant_makeDic(args.vcf)

    # GET CANDIDATE dnSV LIST
    v_parser = md.VariantParser(args.vcf)
    v_parser.parse()
    v_l = v_parser.get_variant_list()
    
    # PHASING CANDIDATE dnSVs
    phased_d = phase_by_germline_info(check_vDic, args.type_info, args.c_bam, args.f_bam, args.m_bam, args.cpu_num, var_count)

    # PRINT PHASED RESULT    
    result_printer(phased_d, args.type_info, args.outname)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = '')
    parser.add_argument('-v', '--vcf', type = str, help = '')
    parser.add_argument('-t', '--type_info', type = str, help = '')
    parser.add_argument('-cb', '--c_bam', type = str, help = '')
    parser.add_argument('-fb', '--f_bam', type = str, help = '')
    parser.add_argument('-mb', '--m_bam', type = str, help = '')
    parser.add_argument('-c', '--cpu_num', type = int, help = '')
    parser.add_argument('-o', '--outname', type = str, help = '')
    args = parser.parse_args()

    import re
    import dnsv_module as md
    import kmer_module as kmer_md

    main(args)
