#!/usr/bin/env python3
import sys
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# read index
def read_index(fr):
    data = {}
    for l in fr:
        ws = l.rstrip('\n').split('\t')
        if len(ws) >= 2:
            data[ws[0]] = int(ws[1])
    return data


# read library list
def read_libs(fr):
    data = {}
    for l in fr:
        ws = l.rstrip('\n').split('\t')
        data[ws[0]] = ws[2]
    return data


# read bed file
def read_bed(fr, use_frequency):
    data = {}
    for l in fr:
        ws = l.rstrip('\n').split('\t')
        if len(ws) < 6:
            continue
        if ws[0] not in data:
            data[ws[0]] = {'+':{}, '-':{}}
        if int(ws[1]) not in data[ws[0]][ws[5]]:
            data[ws[0]][ws[5]][int(ws[1])] = 0
        if use_frequency:
            data[ws[0]][ws[5]][int(ws[1])] += float(ws[3])
        else:
            data[ws[0]][ws[5]][int(ws[1])] += 1
    return data


# get captured fragments:
def get_captured_region(missing, index):
    data = {k:{'+':[], '-':[]} for k in index.keys()}
    for l in missing:
        ws = l.rstrip('\n').split('\t')
        # skip illegal lines
        if len(ws) < 6:
            continue
        # first fragment of each chromosome
        if int(ws[1]) == 0:
            data[ws[0]][ws[5]].append([int(ws[2])])
            continue
        if len(data[ws[0]][ws[5]]) == 0:
            data[ws[0]][ws[5]].append([0])
        data[ws[0]][ws[5]][-1].append(int(ws[1]))
        data[ws[0]][ws[5]].append([int(ws[2])])
    for k in data.keys():
        for s in ['+', '-']:
            if len(data[k][s]) > 0:
                if data[k][s][-1][0] != index[k]:
                    data[k][s][-1].append(index[k])
                else:
                    data[k][s].pop()
    return data


# check if the rNMPs are incorporated in the captured region
def check_rNMPs(rNMPs, region):
    for c in rNMPs.keys():
        for s in ['+', '-']:
            i = 0
            j = 0
            rs = list(sorted(rNMPs[c][s].keys()))
            while i<len(rs) and j<len(region[c][s]):
                if region[c][s][j][1] <= rs[i]:
                    j += 1
                    continue
                elif region[c][s][j][0] > rs[i]:
                    rNMPs[c][s][rs[i]] = [-rNMPs[c][s][rs[i]], False]
                    i += 1
                else:
                    rNMPs[c][s][rs[i]] = [rNMPs[c][s][rs[i]], True]
                    i += 1
            while i < len(rs):
                rNMPs[c][s][rs[i]] = [-rNMPs[c][s][rs[i]], False]
                i += 1
    return rNMPs


# calculate ratios of count and position
def calc_ratio(rNMPs):
    total_counts = 0
    counts_in_region = 0
    positions_in_region = 0
    for s in ['+', '-']:
        for v in rNMPs[s].values():
            if v[1]:
                counts_in_region += abs(v[0])
                positions_in_region += 1
            total_counts += abs(v[0])
    total_positions = len(rNMPs['+']) + len(rNMPs['-'])
    return counts_in_region, total_counts, positions_in_region, total_positions


# calculate rNMP incorporaiton proportion in and out capturable region
def generate_data(libs, bed_folder, missing_folder, index, use_frequency):
    # build list
    groups = {}
    for k, v in libs.items():
        if v not in groups:
            groups[v] = []
        groups[v].append(k)
    # load data
    data = []
    print('loading data...')
    for res, v in groups.items():
        missing_name = f'{missing_folder}/{res}.missings'
        with open(missing_name) as fr:
            captured = get_captured_region(fr, index)
        for lib in v:
            bed_name = f'{bed_folder}/{lib}.bed'
            with open(bed_name) as fr:
                rNMPs = read_bed(fr, use_frequency)
            rNMPs = check_rNMPs(rNMPs, captured)
            for c, rNMPs_in_chrom in rNMPs.items():
                data.append([lib, res, c] + list(calc_ratio(rNMPs_in_chrom)))
            print(f'{lib} done!')
    # convert to df
    df = pd.DataFrame(data, columns=['Library', 'RE', 'Chromosome','Count_in_region',\
         'Total_count', 'Position_in_region', 'Total_position'])
    df['Count_ratio'] = df['Count_in_region']/df['Total_count']
    df['Position_ratio'] = df['Position_in_region']/df['Total_position']
    df = df[['Library', 'RE', 'Chromosome','Count_in_region','Total_count','Count_ratio',\
        'Position_in_region', 'Total_position', 'Position_ratio']]
    print('Data loaded!')
    return df

def main():
    # input command line arguments
    parser = argparse.ArgumentParser(description='Draw bar chart for average coverage in rNMP count and location')
    parser.add_argument('list', type=argparse.FileType('r'), help='List of calculations to be done')
    parser.add_argument('bed_folder', help='Folder of BED file for rNMP incorporation')
    parser.add_argument('missing_folder', help='Folder of BED file for missing parts in reference genome, all missing files with .missing extension should be in this folder')
    parser.add_argument('index', type=argparse.FileType('r'), help='Fasta index file generated by samtools')
    parser.add_argument('-o', default = sys.stdout, type=argparse.FileType('w'), help = 'Output to file')
    parser.add_argument('-f', action='store_true', help='Use the fourth column of the bed file as rNMP incorporation frequency')
    parser.add_argument('--mito_name', default = 'chrM', help = 'Name of mitochondria in reference genome')
    args = parser.parse_args()

    # load data
    index = read_index(args.index)
    libs = read_libs(args.list)

    # calculate ratio of inside/outside region
    df = generate_data(libs, args.bed_folder, args.missing_folder, index, args.f)
    df.to_csv(args.o, sep='\t', index=False)

    print('Done!')



if __name__=='__main__':
    main()
