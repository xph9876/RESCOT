#!/usr/bin/env python3

from resimulation import *
import numpy as np
import argparse
import copy


# class to get cut for optimization
class OptCut(Cut):

    def __init__(self, relist, fasta, min_n, max_n, start, log):
        self.res_c = []
        # the one to remove
        self.cache_res = []
        # the on to add
        self.cache_res_c = []
        self.best_res = []
        self.min_n = min_n
        self.max_n = max_n
        # initialize by random
        if start == '':
            res_idx = np.random.choice(len(relist.res), (min_n + max_n)//2, replace=False)
        # input start
        else:
            res_idx = [int(i) for i in start.split(',')]
        # use get initialize res and res_c
        self.res = [relist.res[i] for i in res_idx]
        for i in relist.res:
            if i not in self.res:
                self.res_c.append(i)
        relist.res = self.res
        super().__init__(fasta, relist)
        self.log = log


    def mutant(self):
        # return :mutant id : 0 for del, 1 for add, 2 for replace
        # re 1 - to del, re 2 - to add
        mutant_id  = 0
        re1 = ''
        re2 = ''
        # when reach the extreme
        if(len(self.res) <= self.min_n):
            # can only add or replace
            if np.random.randint(2) == 0:
                # add
                re1, re2 = self._add_enzyme()
                mutant_id = 1
            else:
                # replace
                re1, re2 = self._replace_enzyme()
                mutant_id = 2
        elif (len(self.res) >= self.max_n):
            # can only delete or replace
            if np.random.randint(2) == 0:
                # delete
                re1, re2 = self._remove_enzyme()
                mutant_id = 0
            else:
                re1, re2 = self._replace_enzyme()
                mutant_id = 2
        else:
            mutant_id = np.random.randint(3)
            if mutant_id == 0:
                # remove
                re1, re2 = self._remove_enzyme()
            elif mutant_id == 1:
                # add
                re1, re2 = self._add_enzyme()
            else:
                # replace
                re1, re2 = self._replace_enzyme()
        self.res.sort()
        return mutant_id,re1,re2
    

    # add an enzyme for mutant
    def _add_enzyme(self):
        i = np.random.randint(len(self.res_c))
        self.cache_res_c = self.res_c[i]
        self.cache_res = ()
        self.res_c.pop(i)
        re2 = self.cache_res_c[0]
        return '', re2


    # remove an enzyme for mutantion
    def _remove_enzyme(self):
        i = np.random.randint(len(self.res))
        self.cache_res = self.res[i]
        self.cache_res_c = ()
        self.res.pop(i)
        re1 = self.cache_res[0]
        return re1, ''


    # replace an enzyme for mutation
    def _replace_enzyme(self):
        # remove
        i = np.random.randint(len(self.res))
        self.cache_res = self.res[i]
        self.res.pop(i)
        re1 = self.cache_res[0]
        # add
        i = np.random.randint(len(self.res_c))
        self.cache_res_c = self.res_c[i]
        self.res_c.pop(i)
        re2 = self.cache_res_c[0]
        return re1, re2


    # accept when reach the threshold of tunneling
    def accept(self):
        if self.cache_res:
            self.res_c.append(self.cache_res[:])
        if self.cache_res_c:
            self.res.append(self.cache_res_c[:])
        # log
        self.log.write('\t' + ','.join([i[0] for i in sorted(self.res)]) + '\taccept\n')


    def decline(self):
        # delete
        if self.cache_res:
            if not self.cache_res_c:
                self.log.write('\t' + ','.join([i[0] for i in sorted(self.res)]) + '\tdecline\n')
            # replace
            else:
                self.log.write('\t' + ','.join([i[0] for i in sorted(self.res + [self.cache_res_c])]) + '\tdecline\n')
                self.res_c.append(self.cache_res_c)
            # add res back
            self.res.append(self.cache_res)
        # add
        else:
            self.log.write('\t' + ','.join([i[0] for i in sorted(self.res + [self.cache_res_c])]) + '\taccept\n')
            self.res_c.append(self.cache_res_c)


    # update best only if have higher coverage
    def update_best(self):
        self.best_res = self.res.copy()


    # generate new cuts when find a new re
    def new_cut(self):
        cuts, _ = self._cut(self.fasta, [self.cache_res_c])
        return cuts


    def export_res(self,filename):
        filename += '.res'
        with open(filename, 'w') as fw:
            for l in sorted(self.best_res):
                fw.write('\t'.join([str(i) for i in l])+'\n')
        print('Optimized restrict enzyme set is out put to {}'.format(filename))


class OptCoverage(Coverage):


    def __init__(self, fragments, min_c, max_c, circulars, log):
        super().__init__(fragments, min_c, max_c, circulars, calc_missings=False)
        self.fragments = fragments.fragments
        self.__best_coverage = self.coverage_overall
        self.log = log
        # likelihood to accept a bad mutantion
        self.beta = 70


    # class to optimize coverage
    def remove_cut(self, enzyme):
        fragments = self.__remove_cut(self.fragments, enzyme)
        # determine whether accept the change or not
        accept, update_best = self.__update(fragments)
        return accept, update_best


    # remove cuts generate by one re from the fragments
    def __remove_cut(self, fragments, enzyme):
        frag_outs = {}
        for chrom, cuts in fragments.items():
            frag_outs[chrom] = [x for x in cuts if x[1]!=enzyme]
        return frag_outs


    def __update(self,fragments):
        accept = False

        # recalculate coverage and overwrite frag_outs
        overall, coverages = self.__calc_overall(fragments)
        self.log.write(str(overall))

        # update best
        if overall > self.__best_coverage:
            update_best = True
            self.__best_coverage = overall
        else:
            update_best = False

        # use stochastic tunneling to jump out the local optimization
        beta = self.beta
        rate_d = (overall - self.coverage_overall)/self.coverage_overall
        accept_prob = np.exp(beta*rate_d)
        if np.random.rand() < accept_prob:
            self.fragments = copy.deepcopy(fragments)
            self.coverages = copy.deepcopy(coverages)
            self.coverage_overall = overall
            accept = True
        return accept, update_best


    # calculate overall coverage
    def __calc_overall(self, fragments):
        coverages= self._calc_coverage(fragments)
        o = self._calculate_overall_coverage(coverages, self.chrom_len)
        return o, coverages


    # to add a set of cuts
    def add_cut(self, new_cut):
        fragments = self.__add_cut(self.fragments, new_cut)
        accept, update_best = self.__update(fragments)
        return accept, update_best


    # add a series of re cuts
    def __add_cut(self, fragments, new_cut):
        frag_outs = {}
        for chrom, cuts in fragments.items():
            frag_outs[chrom] = fragments[chrom] + new_cut[chrom]
            frag_outs[chrom].sort(key=lambda x: x[0])
        return frag_outs
        

    def replace(self, re, cuts):
        fragments1 = self.__remove_cut(self.fragments, re)
        fragments = self.__add_cut(fragments1, cuts)
        accept, update_best = self.__update(fragments)
        return accept, update_best


def main():
    # argument parser
    parser = argparse.ArgumentParser(description='Generate optimized RE set.')
    parser.add_argument('REpool', help='RE candidate pool')
    parser.add_argument('fasta', help='FASTA sequence of genome')
    parser.add_argument('-o', help = 'output file base name, default = FASTA filename')
    parser.add_argument('-l', metavar = 'l', default = 81, type = int, help = 'Minimum suitable genome sequence length, default = 81')
    parser.add_argument('-L', default = 481, type = int, help = 'Maximum suitable genome sequence length, default = 481')
    parser.add_argument('-i', default = 100, type = int, help = 'Number of iterations, default = 100')
    parser.add_argument('-n', default = 1, type = int, help = 'Minimum number of RE in a set, default = 1')
    parser.add_argument('-N', default = 4, type = int, help = 'Maximum number of RE in a set, default = 4')
    parser.add_argument('--circular', nargs='+', default=['chrM'], help='Circular chromosome/plasmids, default=chrM')
    parser.add_argument('--start', default='', type=str, metavar='N,N,N......', help='Index of REs in candidate pool to start optimization with, zero-started index should separated by comma, example:(0,1,3,5)')
    args = parser.parse_args()

    restrict_list_file = args.REpool
    fasta = args.fasta
    if not args.o:
        args.o = fasta
    max_capture = args.L
    min_capture = args.l
    iterations = args.i
    max_number = args.N
    min_number = args.n

    # initialize
    # open log
    logw = open(args.o + '.log', 'w')
    re_list= REList(restrict_list_file)
    print(f'Totally {len(re_list.res)} restrict enzymes are found, start optimization!')
    optcut = OptCut(re_list,fasta, min_n=args.n, max_n=args.N, start=args.start, log=logw)
    fragments = Fragments(optcut)
    optcoverage = OptCoverage(fragments, min_capture, max_capture, args.circular, log=logw)

    # write first coverage to log
    logw.write(f'{optcoverage.coverage_overall}\t' + ','.join([x[0] for x in optcut.res]) + '\taccept\n')

    # optimize
    for i in range(iterations):
        s = ''
        s += 'Iteration {} -- '.format(i+1)
        mutant_id, re_del, re_add = optcut.mutant()
        # delete
        if mutant_id == 0:
            s += 'delete {} --'.format(re_del)
            accept, update_best = optcoverage.remove_cut(re_del)
            if accept:
                optcut.accept()
                if update_best:
                    optcut.update_best()
                s += 'accept!'
            else:
                optcut.decline()
                s += 'decline!'
        # add
        elif mutant_id == 1:
            s += 'add {} -- '.format(re_add)
            cuts = optcut.new_cut()
            accept, update_best = optcoverage.add_cut(cuts)
            if accept:
                optcut.accept()
                if update_best:
                    optcut.update_best()
                s += 'accept!'
            else:
                optcut.decline()
                s += 'decline!'
        # replace
        else:
            s += 'replace {} to {} -- '.format(re_del,re_add)
            cuts = optcut.new_cut()
            accept, update_best = optcoverage.replace(re_del,cuts)
            if accept:
                optcut.accept()
                if update_best:
                    optcut.update_best()
                s += 'accept!'
            else:
                optcut.decline()
                s += 'decline!'
        print(s)

    # output
    # res
    optcut.export_res(args.o)
    res = REList(args.o + '.res')

    # cut
    cut = Cut(fasta, res)
    cut.export(args.o)

    # fragment
    fragments = Fragments(cut)
    del cut
    fragments.export(args.o)

    # coverage
    coverage = Coverage(fragments, min_capture, max_capture, args.circular)
    del fragments
    coverage.export(args.o)



if __name__ == '__main__':
    main()
