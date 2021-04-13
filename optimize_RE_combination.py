#!/usr/bin/env python3
import argparse
import resimulation as res
import sys
import numpy as np
from collections import defaultdict
from tqdm import tqdm


# class for random restriction enzyme set
class RandomRESet(res.REList):
    def rand_samp(self, num):
        idx = np.random.choice(np.arange(len(self.res)), num, replace=False)
        idx.sort()
        self.idx = tuple(idx)
        return tuple(idx)
    def output(self, output):
        with open(output, 'w') as fw:
            for i in self.idx:
                fw.write('\t'.join([str(k) for k in self.res[i]]) + '\n')


def main():
    # argument parser
    parser = argparse.ArgumentParser(description='Generate optimized RE combination')
    parser.add_argument('REpool',help='RE candidate pool')
    parser.add_argument('FASTA', help='Sequence file in FASTA format')
    parser.add_argument('-o', help = 'output file base name, default = FASTA filename')
    parser.add_argument('-l', metavar = 'l', default = 81, type = int, help = 'Minimum suitable genome sequence length, default = 81')
    parser.add_argument('-L', default = 481, type = int, help = 'Maximum suitable genome sequence length, default = 481')
    parser.add_argument('-n', default = 1, type = int, help = 'Minimum number of RE in a set, default = 1')
    parser.add_argument('-N', default = 4, type = int, help = 'Maximum number of RE in a set, default = 4')
    parser.add_argument('-m', type=int, default=2, help='Number of restriction enzyme sets in the combination, default=2')
    parser.add_argument('-i', default = 100, type = int, help = 'Number of random RE sets generated, default = 100')
    parser.add_argument('-c', type=int, default=1000, help='Number of iterations for combination, default=200000')
    parser.add_argument('--circular', nargs='+', default=['chrM'], help='Circular chromosome/plasmids, default=chrM')
    parser.add_argument('--combination_only', action='store_true', help='Only do combiantion')
    args = parser.parse_args()

    # output file basename:
    if args.o == '':
        OUTPUT_FILE_BASE = args.FASTA
    else:
        OUTPUT_FILE_BASE = args.o

    # read restriction enzyme
    all_re = RandomRESet(args.REpool)

    # build the RE sets
    resets = {}
    if not args.combination_only:
        # random sampling for restriction enzyme set
        print('Start to generate restriction enzyme cut!')
        hist = {}
        for i in tqdm(range(args.i), desc='Generating RE set'):
            num = np.random.randint(args.n, args.N+1)
            sampled = all_re.rand_samp(num)
            # append to reset with out replication
            if sampled not in hist:
                hist[sampled] = 1
                # output
                output = OUTPUT_FILE_BASE + '_{:d}.res'.format(i)
                all_re.output(output)
                relist = res.REList(output)
                reset = tuple(sorted([x[0] for x in relist.res]))
                if reset not in list(resets.values()):
                    resets[i] = reset
                    cut = res.Cut(args.FASTA, relist)
                    fragments = res.Fragments(cut)
                    coverage = res.Coverage(fragments, min_c=args.l, max_c=args.L, circulars=args.circular, calc_missings=True)
                    coverage.export_missing(f'{OUTPUT_FILE_BASE}_{i}.missings', silent=True)
    else:
        # reading from generated file
        print('Reading RE sets...')
        for i in range(args.i):
            output = OUTPUT_FILE_BASE + '_{:d}.res'.format(i)
            try:
                relist = res.REList(output)
            except FileNotFoundError:
                continue
            reset = tuple(sorted([x[0] for x in relist.res]))
            if reset not in list(resets.values()):
                resets[i] = reset

    # random sampling for combination
    least_missing = 2**32
    best_comb = None
    print('Start random combination!')
    with open(args.o + '.log', 'w') as logw:
        logw.write('\t'.join(['res{:d}'.format(x+1) for x in range(args.m)]) + '\tMissing\n')
        hist = {}
        for it in tqdm(range(args.c), desc='Combination'):
            # random choose n sets
            idx = np.random.choice(list(resets.keys()), args.m, replace=False)
            idx = tuple(sorted(idx))
            # skip duplicates
            if not idx in hist:
                hist[idx] = 1
            else:
                continue
            # skip duplicates
            dup = []
            isdup = False
            for i in idx:
                if not sorted(resets[i]) in dup:
                    dup.append(sorted(resets[i]))
                else:
                    isdup=True
                    break
            if isdup:
                continue

            # get intersect
            bed_dict = {'+' : defaultdict(lambda :defaultdict(int)), '-':defaultdict(lambda : defaultdict(int))}
            for i in idx:
                with open(OUTPUT_FILE_BASE + '_{}.missings'.format(i)) as fr:
                    for l in fr:
                        ws = l.rstrip('\n').split('\t')
                        if len(ws) < 6:
                            continue
                        bed_dict[ws[5]][ws[0]][int(ws[1])] += 1
                        bed_dict[ws[5]][ws[0]][int(ws[2])] -= 1

            # calculate length of missing parts of combinations
            missing = 0
            in_intersect= False
            for s in ['+', '-']:
                for chrom, v in bed_dict[s].items():
                    if not v:
                        continue
                    count = 0
                    for i in sorted(v.keys()):
                        count += v[i]
                        assert count <= args.m, f'Error in count missings of combination {idx} in chromosome {chrom}({s})!'
                        if in_intersect:
                            if count < args.m:
                                missing += i
                                in_intersect = False
                        else:
                            if count == args.m:
                                missing -= i
                                in_intersect = True
            logw.write('\t'.join([','.join([r for r in resets[x]]) for x in idx]) + '\t{:d}\n'.format(missing))

            if least_missing > missing:
                least_missing = missing
                best_comb = idx[:]

    # print
    print('Best combination is ' + ', '.join([(OUTPUT_FILE_BASE + '_' + str(i)) for i in best_comb]) + \
            '. Which has {} bp missing totally.'.format(least_missing))

if __name__ == '__main__':
    main()

