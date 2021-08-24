#!/home/liujf/WORKSPACE/zhuoy/bin/pypy
'''
Description:
Author: Zhuo Yue
Date: 2021-06-09 15:02:54
LastEditors: Zhuo Yue
LastEditTime: 2021-08-24 16:23:40
Calls:
Called By:
FilePath: \CDpan\bin\ex.py
'''

import time
import copy
import re
import os
import sys


# Screening threshold control function
def Rate(pairwise_list, threshold):
    # threshold = 0.95
    if float(pairwise_list[9]) / float(pairwise_list[6]) >= threshold:
        # if float(pairwise_list[9]) / float(pairwise_list[1]) >= threshold:
        return True
    return False


def Splicing(pairwise_list_list, threshold):
    # When there are multiple matches, try to merge
    # The sum of matched bases is less than the threshold, no need to merge
    pairwise_new = copy.deepcopy(pairwise_list_list[0])
    pairwise_new[9] = 0
    for line in pairwise_list_list:
        pairwise_new[9] += int(line[9])
    if not Rate(pairwise_new, threshold):
        return False

    # Bubble sort, adjust the order according to the matching position
    # Quick sort is troublesome and don’t want to write, anyway, it’s not useful to get this part
    counter = True
    while counter:
        counter = False
        for i in range(len(pairwise_list_list) - 1):
            if int(pairwise_list_list[i][7]) > int(pairwise_list_list[i+1][7]):
                pairwise_list_list[i], pairwise_list_list[i +
                                                          1] = pairwise_list_list[i+1], pairwise_list_list[1]
                counter = True

    # Try to merge
    pairwise_new = copy.deepcopy(pairwise_list_list[0])
    have_merged = False
    for i in range(1, len(pairwise_list_list)):
        if int(pairwise_list_list[i][7]) - int(pairwise_new[8]) < 100:
            if int(pairwise_list_list[i][7]) - int(pairwise_new[8]) > 0:
                pairwise_new[8] = pairwise_list_list[i][8]
                pairwise_new[9] = str(
                    int(pairwise_new[9]) + int(pairwise_list_list[i][9]))
                have_merged = True
            else:
                pass
        else:
            pairwise_new = copy.deepcopy(pairwise_list_list[i])
    # No need to compare thresholds for unmerged
    if have_merged:
        return Rate(pairwise_new, threshold)
    else:
        return False


def Comparison(target_string, query_dict, threshold):
    # Match
    if target_string not in query_dict:
        return '0'

    for value in query_dict.get(target_string).values():
        if len(value) == 1:
            if Rate(value[0], threshold):
                return '1'
        elif len(value) > 1:
            for value_value in value:
                if Rate(value_value, threshold):
                    return '1'

            if Splicing(value, threshold):
                return '1'

    return '0'


def Arrange(paf_file, seq_list, location_dict, threshold):
    pad = {}
    with open(paf_file) as f:
        for line in f:
            line = line.rstrip('\n').split('\t')[:12]
            if line[5] not in pad:
                pad[line[5]] = {line[0]: [line]}
                continue

            if line[0] not in pad.get(line[5]):
                pad[line[5]][line[0]] = [line]
            else:
                pad[line[5]][line[0]].append(line)

    for i in range(len(seq_list)):
        com = Comparison(seq_list[i][0], pad, threshold)
        status = location_dict.get(seq_list[i][0], '')
        seq_list[i].append(f'{com}{status}')

    return seq_list


# def Aexist(coverage_file, seq_list, location_dict):
#     coverage = {}
#     with open(coverage_file) as f:
#         for line in f:
#             if re.match("^#", line):
#                 continue

#             line = line.rstrip('\n').split('\t')
#             if float(line[5]) > 95:
#                 coverage[line[0]] = 1

#     for i in range(len(seq_list)):
#         status = location_dict.get(seq_list[i][0], '')
#         com = coverage.get(seq_list[i][0], 0)
#         seq_list[i].append(f'{com}{status}')

#     return seq_list


def Location(file_path_paf_string, file_path_location_string):
    breed_individual = re.search(
        '[^/]+\/[^/]+\.all\.paf', file_path_paf_string).group().replace('.all.paf', '')
    file_location = f'{file_path_location_string}{breed_individual}'
    # TODO 1/0 check?
    location = {}

    with open(f'{file_location}/2a.name') as f:
        for line in f:
            seq = line.rstrip('\n').split(' ')

            status = 'A'
            contig = seq[0]
            chr = seq[2]
            left = (int(seq[4]) + int(seq[3])) // 2
            right = (int(seq[6]) + int(seq[5])) // 2

            if (left > right):
                left, right = right, left

            if (abs(right - left) > 500):
                status = 'N'
                chr = ''
                left = ''
                right = ''

            location[contig] = f':{status}:{chr}:{left}:{right}'

    with open(f'{file_location}/2bleft.name') as f:
        for line in f:
            seq = line.rstrip('\n').split(' ')

            status = 'L'
            contig = seq[0]
            chr = seq[2]
            left = (int(seq[4]) + int(seq[3])) // 2
            right = ''

            location[contig] = f':{status}:{chr}:{left}:{right}'

    with open(f'{file_location}/2bright.name') as f:
        for line in f:
            seq = line.rstrip('\n').split(' ')

            status = 'R'
            contig = seq[0]
            chr = seq[2]
            left = ''
            right = (int(seq[4]) + int(seq[3])) // 2

            location[contig] = f':{status}:{chr}:{left}:{right}'

    with open(f'{file_location}/3.name') as f:
        for line in f:
            seq = line.rstrip('\n').split(' ')

            status = 'U'
            contig = seq[0]
            chr = seq[1]
            left = ''
            right = ''

            location[contig] = f':{status}:{chr}:{left}:{right}'

    with open(f'{file_location}/4.name') as f:
        for line in f:
            seq = line.rstrip('\n').split(' ')

            status = 'N'
            contig = seq[0]
            chr = ''
            left = ''
            right = ''

            location[contig] = f':{status}:{chr}:{left}:{right}'

    with open(f'{file_location}/5.name') as f:
        for line in f:
            seq = line.rstrip('\n').split(' ')

            status = 'N'
            contig = seq[0]
            chr = ''
            left = ''
            right = ''

            location[contig] = f':{status}:{chr}:{left}:{right}'

    return location


debug = False

for threshold in [x / 100 for x in range(70, 100, 5)]:
    # if debug:
    #     if (threshold > 0.71):
    #         continue
    # Start of the main program
    start_time = time.time()
    ''''
    if len(sys.argv) != 5:
        print('Error: Number of wrong parameters.')
        sys.exit(1)
    fasta, arrange, location, output = sys.argv[1:5]
    '''
    fasta = '/home/liujf/WORKSPACE/zhuoy/test/dispensable_genome.fasta.length'

    # arrange = '/home/liujf/WORKSPACE/zhuoy/test/arrange/'
    # arrange = '/home/liujf/WORKSPACE/zhuoy/test/test/'

    location = '/home/liujf/WORKSPACE/zhuoy/test/location/'

    output = "/home/liujf/WORKSPACE/zhuoy/test/compare"
    if debug:
        output = '/home/liujf/WORKSPACE/zhuoy/test/test'

    # aexist = '/home/liujf/WORKSPACE/zhuoy/test/aexist/'
    # aexist = '/home/liujf/WORKSPACE/zhuoy/test/test/'

    minimap = '/home/liujf/WORKSPACE/zhuoy/test/minimap/'
    if debug:
        minimap = '/home/liujf/WORKSPACE/zhuoy/test/test/'

    seq = []
    with open(fasta) as f:
        for line in f:
            seq.append([line.rstrip('\n').split('\t')[0]])

    header = ['Contig']
    os.chdir(minimap)
    file_list = os.popen('find')
    for file_string in file_list:
        file = re.search('[^/]+\.all\.paf', file_string)
        if file:
            file = file.group()
            # TODO 先做一步contigs存在性的判定
            loc = Location(file_string.rstrip('\n'), location)
            seq = Arrange(file_string.rstrip('\n'), seq, loc, threshold)
            # seq = Aexist(file_string.rstrip('\n'), seq, loc)
            header.append(file.replace('.all.paf', ''))

    f = open(f"{output}_{threshold}.txt", 'w+')
    f.write(' '.join(header) + '\n')
    for line in seq:
        f.write(' '.join(line) + '\n')
    f.close()

    end_time = time.time()
    print(f"It took {end_time-start_time:.2f} seconds to compute.\n")
