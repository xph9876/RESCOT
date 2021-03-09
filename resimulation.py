#!/usr/bin/env python3
import sys
import re
import argparse
import numpy as np
from collections import defaultdict


class REList(object):

    def __init__(self, re_file):
        self._DEGENRATE_BASE={'R':'[AG]', 'Y':'[CT]','M':'[AC]', 'K':'[GT]', 'S':'[CG]', 'W':'[AT]', 'H':'[ACT]', 'B':'[CGT]', 'D':'[AGT]', 'V':'[ACG]', 'N':'[ACGT]'}
        self.res = self._get_res(re_file, self._DEGENRATE_BASE)

    def _get_res(self, re_file, DEGENRATE_BASE):
        res = []
        with open(re_file) as fr:
            for l in fr:
                words = l.rstrip('\n').split('\t')
                assert len(words) == 3, 'Wrong format restrict enzyme list file, the right one is 3 columns separated with TAB, such as:\nEcoRV\tGATATC\t3'
                # remove spaces
                words[1] = words[1].replace(' ','')
                temp = ''
                for c in words[1].upper():
                    if c in DEGENRATE_BASE.keys():
                        temp += DEGENRATE_BASE[c]
                    else:
                        temp += c
                res.append((words[0], temp, int(words[2])))
        return res


class Cut(object):

    # class for restriction enzyme cut
    fasta = ''
    res = []
    # cuts[chrom] = [[pos, re]....]
    chrom_len = {}


    def __init__(self, fasta, relist):
        self.res = relist.res
        self.fasta = fasta
        self.cuts, self.chrom_len = self._cut(fasta, self.res)


    def _cut(self,fasta, res):
        # search the cut site by chromosome
        coordinate_list = defaultdict(list)
        chromosome_length = {}
        chromosome = ''
        seq = ''
        with open(fasta) as fr:
            for line in fr:
                line = line.rstrip('\n')
                if line == '':
                    continue
                elif line[0] == '>':
                    #add information of chromosome and search cut sites when meet the ends of chromosome
                    for i in range(len(res)):
                        # coordinate_list = [cut sets]
                        # position of last base of first half (start with 1)
                        coordinate_list[chromosome] += [[m.start()+ int(res[i][2]), res[i][0]] for m in re.finditer(res[i][1],seq)]
                    chromosome_length[chromosome] = len(seq)
                    #clear information
                    chromosome = line[1:]
                    seq = ''
                else:
                    seq += line.upper()
        #search cut site for last chromosome
        for i in range(len(res)):
            coordinate_list[chromosome] += [[m.start()+ int(res[i][2]), res[i][0]] for m in re.finditer(res[i][1],seq)]
        chromosome_length[chromosome] = len(seq)
        # remove empty
        if '' in chromosome_length:
            del chromosome_length['']
        if '' in coordinate_list:
            del coordinate_list['']
        assert len(chromosome_length)>0 , 'Cannot find any sequences in fasta file'
        # sort
        for i in coordinate_list.keys():
            coordinate_list[i] = sorted(coordinate_list[i], key=lambda x:x[0])
        return coordinate_list, chromosome_length


    # output cuts to file
    def export(self, filename):
        filename += '.cuts'
        # chromosome     position    restrict_enzyme
        with open(filename,'w') as w:
            for chrom in self.chrom_len.keys():
                for cut in self.cuts[chrom]:
                    w.write(f'{chrom}\t{cut[0]}\t{cut[1]}\n')
        print (f'Output cut site file: {filename}!')


class Fragments(object):

    def __init__(self, cuts):
        self.res = cuts.res
        self.chrom_len = cuts.chrom_len
        self.fragments = self._get_fragments(cuts.cuts)


    def _get_fragments(self, cuts):
        fragments = cuts
        for chrom in cuts.keys():
            fragments[chrom] = [[0, 'start']] + fragments[chrom]
            fragments[chrom].append([self.chrom_len[chrom], 'end'])
        return fragments


    def export(self, filename):
        filename += '.fragments'
        with open(filename, 'w') as fw:
            for chrom, fms in self.fragments.items():
                for i in range(0, len(fms)-1):
                    fw.write(f'{chrom}\t{fms[i][0]}\t{fms[i+1][0]}\t{fms[i][1]}\t{fms[i+1][1]}\t{fms[i+1][0]-fms[i][0]}\n')
        print ('Fragments are output to {}'.format(filename))


class Coverage(Fragments):


    def __init__(self, fragments, min_c, max_c, circulars, calc_missings=False):
        self.min_capture = min_c
        self.max_capture = max_c
        self.chrom_len = fragments.chrom_len
        self.calc_missings = calc_missings
        self.circulars = circulars
        self._calc(fragments.fragments)


    def _calc(self, fragments):
        self.coverages = self._calc_coverage(fragments)
        if self.calc_missings:
            self.missings_fwd, self.missings_rev = self._calc_missings(fragments)
        self.coverage_overall = self._calculate_overall_coverage(self.coverages, self.chrom_len)


    def _calc_coverage(self, fragments):
        # load data
        chromosome_length = self.chrom_len
        min_capture = self.min_capture
        max_capture = self.max_capture
        missings = {k:0 for k in chromosome_length.keys()}
        for k, v in fragments.items():
            for i in range(len(v)-1):
                l = v[i+1][0] - v[i][0]
                # first missing part
                missings[k] += min(min_capture, l)
                # second missing part
                missings[k] += max(0, l-max_capture)
            # correct circulars
            if k in self.circulars:
                # no cut in the circular region, nothing captured
                if len(v) == 2:
                    missings[k] = chromosome_length[k]
                    continue
                # remove first and last
                for i in [1, -1]:
                    l = v[i][0]-v[i-1][0]
                    missings[k] -= min(min_capture, l)
                    missings[k] -= max(0, l-max_capture)
                # add combined
                l = v[1][0] + v[-1][0] - v[-2][0]
                missings[k] += min(min_capture, l)
                missings[k] += max(0, l-max_capture)
        # calc coverage from missings
        coverage = {k:1-missings[k]/chromosome_length[k] for k in missings.keys()}
        return coverage


    def _calc_missings(self, fragments):
        # load data
        chromosome_length = self.chrom_len
        min_capture = self.min_capture
        max_capture = self.max_capture
        missings_fwd = {k:[] for k in chromosome_length.keys()}
        missings_rev = {k:[] for k in chromosome_length.keys()}
        for k, v in fragments.items():
            # correct circulars
            if k in self.circulars:
                # no cut in the circular region, nothing captured
                if len(v) == 2:
                    missings_fwd[k].append([0, chromosome_length[k]])
                    missings_rev[k].append([0, chromosome_length[k]])
                    continue
                # deal with the first fragment
                cache_fwd = []
                cache_rev = []
                self._calc_missings_helper(k, missings_fwd, missings_rev, v[-2][0]-v[-1][0], v[1][0], min_capture, max_capture)
                # move all missing part before 0 to last
                # forward
                while missings_fwd[k] and missings_fwd[k][0][0]<0:
                    if missings_fwd[k][0][1] <= 0:
                        cache_fwd.append([chromosome_length[k]+missings_fwd[k][0][0], chromosome_length[k]+missings_fwd[k][0][1]])
                        missings_fwd[k].pop(0)
                    else:
                        cache_fwd.append([chromosome_length[k]+missings_fwd[k][0][0], chromosome_length[k]])
                        missings_fwd[k][0][0] = 0
                # reverse
                while missings_rev[k] and missings_rev[k][0][0]<0:
                    if missings_rev[k][0][1] <= 0:
                        cache_rev.append([chromosome_length[k]+missings_rev[k][0][0], chromosome_length[k]+missings_rev[k][0][1]])
                        missings_rev[k].pop(0)
                    else:
                        cache_rev.append([chromosome_length[k]+missings_rev[k][0][0], chromosome_length[k]])
                        missings_rev[k][0][0] = 0
                # deal with the rest
                for i in range(1, len(v)-2):
                    self._calc_missings_helper(k, missings_fwd, missings_rev, v[i][0], v[i+1][0], min_capture, max_capture)
                # append the last
                missings_fwd[k] += cache_fwd
                missings_rev[k] += cache_rev
            # no circular
            else:
                for i in range(len(v)-1):
                    self._calc_missings_helper(k, missings_fwd, missings_rev, v[i][0], v[i+1][0], min_capture, max_capture)
        return missings_fwd, missings_rev


    # helper function of _calc_missings()
    def _calc_missings_helper(self, chrom, missings_fwd, missings_rev, start, end, min_capture, max_capture):
        l = end-start
        # forward strand
        # first half on the left
        if min_capture:
            if missings_fwd[chrom] and missings_fwd[chrom][-1][1] == start:
                missings_fwd[chrom][-1][1] += min(min_capture, l)
            else:
                missings_fwd[chrom].append([start, start+min(min_capture, l)])
        # second half on the left
        if l > max_capture:
            missings_fwd[chrom].append([start+max_capture, end])
        # reverse strand
        # first half on the left
        if l > max_capture:
            if missings_rev[chrom] and missings_rev[chrom][-1][1] == start:
                missings_rev[chrom][-1][1] += (l-max_capture)
            else:
                missings_rev[chrom].append([start, start+l-max_capture])
        # second half on the left
        if min_capture:
            missings_rev[chrom].append([end-min(min_capture, l), end])


    def _calculate_overall_coverage(self, chromosome_coverage, chromosome_length):
        cover = 0
        l = 0
        for i in chromosome_coverage.keys():
            cover += chromosome_coverage[i]*chromosome_length[i]
            l += chromosome_length[i]
        return cover/l


    def export(self, basename, silent=False):
        if self.calc_missings:
            self.export_missing(basename +'.missings')
        self.export_coverage(basename + '.coverage')
        if not silent:
            print ('Overall coverage: {:.4f}'.format(self.coverage_overall))


    def export_missing(self,filename, silent=False):
        with open(filename, 'w') as fw:
            for chrom in self.chrom_len.keys():
                for l in self.missings_fwd[chrom]:
                    fw.write('\t'.join([chrom, str(l[0]), str(l[1]), '.', '.', '+']) + '\n')
                for l in self.missings_rev[chrom]:
                    fw.write('\t'.join([chrom, str(l[0]), str(l[1]), '.', '.', '-']) + '\n')
        if not silent:
            print('Missing parts are output to '+ filename)


    def export_coverage(self,filename, silent=False):
        with open(filename, 'w') as fw:
            fw.write(f'OVERALL\t{self.coverage_overall:.4f}\n')
            for k,v in self.coverages.items():
                fw.write(f'{k}\t{v:.4f}\n')
        if not silent:
            print('Coverages for each chromosome are export to '+  filename)



def main():
    # input command line arguments
    parser = argparse.ArgumentParser(description='Calculate the capture percent and export missing part when sequencing the fragments made by restrict enzyme.', prog='python calculate_uncapture_sequence.py')
    parser.add_argument('res', help='list of restrict enzyme cut sites')
    parser.add_argument('fasta', help='FASTA sequence of genome')
    parser.add_argument('-l', metavar = 'l', default = 81, type = int, help = 'length of the former missing part, every fragments contain this part, default = 81')
    parser.add_argument('-nl', action = 'store_true', help = 'no former missing part, default = False')
    parser.add_argument('-L', default = 481, type = int, help = 'max capture length, default = 481')
    parser.add_argument('-o', default = '', help = 'output file base name, default = FASTA filename')
    parser.add_argument('--circular', nargs='+', default=['chrM'], help='Circular chromosome/plasmids, default=chrM')
    args = parser.parse_args()

    # input informations
    if not args.o:
        args.o = args.fasta
    max_capture = args.L
    if args.nl:
        min_capture = 0
    else:
        min_capture = args.l

    # get restriction enzymes
    res = REList(args.res)
    print(f'Found {len(res.res)} restrict enzymes. Start searching!')

    # generate cuts
    cut = Cut(args.fasta, res)
    cut.export(args.o)

    # generate fragments
    fragments = Fragments(cut)
    del cut
    fragments.export(args.o)
    
    # calculate coverage
    coverage = Coverage(fragments, min_capture, max_capture, args.circular, calc_missings=True)
    del fragments
    coverage.export(args.o)

    print('Done!')



if __name__=='__main__':
    main()
