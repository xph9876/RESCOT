#!/usr/bin/env python3
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# read logs
def read_log(frs, genome_size):
    data = []
    for fr in frs:
        # skip header
        fr.readline()
        # initialize
        i = 0
        mn = 2**32
        for l in fr:
            ws = l.rstrip('\n').split('\t')
            if len(ws) < 2:
                continue
            coverage = float(ws[-1])
            mn = min(coverage, mn)
            data.append([i, fr.name.split('.')[0], 1-(mn/genome_size/2)])
            i += 1
    data = pd.DataFrame(data, columns=['Iteration', 'Experiment', 'Coverage'])
    return data


# get genome size
def get_genome_size(fr):
    s = 0
    for l in fr:
        ws = l.rstrip('\n').split('\t')
        if len(ws) < 3:
            continue
        s += int(ws[1])
    return s


def main():
    # input command line arguments
    parser = argparse.ArgumentParser(description='Draw line chart for bestn log')
    parser.add_argument('log', type=argparse.FileType('r'), nargs='+', help='log files of best n, each file forms a line')
    parser.add_argument('index', type=argparse.FileType('r'), help='Fasta index file of reference genome.')
    parser.add_argument('-o', default = 'bestn_analysis', help = 'Output basename, default = bestn_analysis')
    args = parser.parse_args()

    # load data
    genome_size = get_genome_size(args.index)
    df = read_log(args.log, genome_size)

    # draw bar plots
    sns.set(style='ticks', font_scale=2)
    fig, ax = plt.subplots(figsize=(6,6))
    plt.subplots_adjust(left=0.18, right=0.95, top=0.95, bottom=0.14)
    sns.lineplot(x='Iteration', y='Coverage', hue='Experiment', data=df, palette='Set2', linewidth=3, ax=ax)
    # ax.get_legend().remove()
    sns.despine()
    # x and y limit
    plt.xlim([0, df.Iteration.max()])
    # label
    plt.ylabel(f'Coverage' )
    plt.xlabel('Iteration')
    plt.savefig(f'{args.o}_line.png')

    print('Done!')



if __name__=='__main__':
    main()
