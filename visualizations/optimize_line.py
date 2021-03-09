#!/usr/bin/env python3
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# read index
def read_log(fr):
    data = []
    i = 0
    last_accepted = 0
    mx = 0
    for l in fr:
        ws = l.rstrip('\n').split('\t')
        if len(ws) < 3:
            continue
        coverage = float(ws[0])
        data.append([i, 'Calculated', coverage])
        # accepted
        if ws[2] == 'accept':
            last_accepted = coverage
        data.append([i, 'Accepted', last_accepted])
        # maximum
        mx = max(mx, coverage)
        data.append([i, 'Maximum', mx])
        i += 1
    data = pd.DataFrame(data, columns=['Iteration', 'Line', 'Coverage'])
    return data


def main():
    # input command line arguments
    parser = argparse.ArgumentParser(description='Draw line chart for optimization log')
    parser.add_argument('tsv', type=argparse.FileType('r'), help='TSV file of optimization log')
    parser.add_argument('-o', default = 'optimization_analysis', help = 'Output basename, default = capturable_analysis')
    args = parser.parse_args()

    # load data
    df = read_log(args.tsv)

    # draw bar plots
    sns.set(style='ticks', font_scale=2)
    fig, ax = plt.subplots(figsize=(6,6))
    plt.subplots_adjust(left=0.18, right=0.95, top=0.95, bottom=0.14)
    sns.lineplot(x='Iteration', y='Coverage', hue='Line', hue_order=['Calculated','Accepted', 'Maximum'], data=df, ax=ax)
    ax.get_legend().remove()
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
