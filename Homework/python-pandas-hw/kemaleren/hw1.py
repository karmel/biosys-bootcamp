"""Bootcamp homework

Author: Kemal Eren

"""

from __future__ import division

import numpy as np
import itertools
import pandas as pd


class ExpressionAnalyzer(object):
    """Imports and analyzes gene expression data."""
    def __init__(self):
        self.answers_file = open('./answers.txt', 'w')
        self.answers_file.write('Bootcamp homework\n')
        self.answers_file.write('Kemal Eren\n\n')

    def _cell_lines(self):
        return set([k.split('_')[0] for k in self.df.keys()
                    if 'hrs' in k])

    def _write_answer(self, letter, string):
        self.answers_file.write('Question {}\n'.format(letter))
        self.answers_file.write('----------\n')
        self.answers_file.write(string)
        self.answers_file.write('\n\n')

    def import_file(self, filename):
        """Import data from file.

        Expects data with named samples as columns and genes as rows.

        """
        df = pd.read_csv(filename, sep="\t")

        # drop 'call' columns
        subset_keys = list(k for k in df.keys() if k[:4] != 'call')
        self.df = df[subset_keys]

    def do_a(self):
        """Question A: get number of unique genes.

        Assume X is the same as X_i, X_r, etc.
        """
        n_unique = len(set(n.split('_')[0]
                           for n in self.df['Gene Accession Number']))
        answer1 = ('There are {} unique genes according to'
                   ' accession number.'.format(n_unique))
        answer2 = ('There are {} unique genes according to'
                   ' gene description.'.format(
                       len(set(self.df['Gene Description']))))

        self._write_answer('A', '\n'.join([answer1, answer2]))

    def do_b(self):
        within_cell_line_corrs = {}
        for cl in self._cell_lines():
            cl_keys = list(k for k in self.df.keys() if cl in k)
            pairs = list(itertools.combinations(cl_keys, 2))
            corrs = dict(((a.split('_')[1], b.split('_')[1]),
                          abs(self.df[a].corr(self.df[b])))
                         for a, b in pairs)
            maxkey = max(corrs, key=corrs.get)
            within_cell_line_corrs[cl] = maxkey
        answer = []
        for key, val in within_cell_line_corrs.items():
            h1, h2 = val
            answer.append('Cell line {} is most correlated'
                          ' in hours {} and {}.'.format(
                              key, h1, h2))

        self._write_answer('B', '\n'.join(answer))

    def do_c(self):
        """Cell line similarity.

        Assumes we are interested in base expression similarity, not
        reaction to treatment.

        """
        cl_keys = list(k + "_0_hrs" for k in self._cell_lines())
        pairs = list(itertools.combinations(cl_keys, 2))
        corrs = dict(((a.split('_')[0], b.split('_')[0]),
                      abs(self.df[a].corr(self.df[b])))
                     for a, b in pairs)
        line1, line2 = max(corrs, key=corrs.get)
        answer = ('Cell lines {} and {} are the most'
                  ' similar.'.format(line1, line2))
        self._write_answer('C', answer)

    def do_d(self):
        """unchanging genes: lowest variance across all cell lines"""
        n_take = 10
        idxs = self.df.var(axis=1).argsort()[:n_take]
        most_unchanging_genes = list(self.df["Gene Accession Number"][idxs])
        answer = 'The {} lowest variance genes in all lines:\n{}'.format(
            n_take, '\n'.join(most_unchanging_genes))
        self._write_answer('D', answer)

    def _fold_change(self, a, b):
        """Returns fold change from a to b.

        If any values are negative, assumes they are 0 (i.e. no
        expression). If a is 0 and b is nonzero, returns a large
        number. If a is nonzero and b is zero, returns a large
        negative number.

        Converts fractional change to negative inverse.

        """
        a = self.df[a].copy()
        b = self.df[b].copy()
        a[a < 0] = 0
        b[b < 0] = 0
        vals = b / a

        idxs = vals < 1
        vals[idxs] = (-1 / vals[idxs])

        vals[(a == 0) & (b == 0)] = 1
        vals[(a == 0) & (b > 0)] = 1000
        vals[(a > 0) & (b == 0)] = -1000

        return vals

    def do_e(self):
        """two-fold higher expression at 24 hours for all for cell
        types"""
        bools = list(self._fold_change(cl + '_0_hrs',
                                       cl + '_24_hrs').abs() >= 2
                     for cl in self._cell_lines())
        bools = reduce(np.logical_and, bools)
        higher_expression = self.df["Gene Accession Number"][bools]
        answer = ('There are {} genes with two-fold higher expression at'
                  ' 24 hours in all cell lines: {}'.format(
                      len(higher_expression), ', '.join(higher_expression)))
        self._write_answer('E', answer)

    def do_f(self):
        """differentially regulated genes"""
        bools = self._fold_change('U937_0_hrs', 'HL60_0_hrs').abs() >= 2
        diff_regulated = self.df["Gene Accession Number"][bools]
        answer = ('There are {} genes differentially regulated in'
                  ' HL60 vs. U937 at 0 hours:'
                  ' {}'.format(len(diff_regulated), ', '.join(diff_regulated)))
        self._write_answer('F', answer)
        diff_regulated = set(g.split('_')[0] for g in diff_regulated)
        with open('./diff_regulated.txt', 'w') as f:
            f.write('\n'.join(diff_regulated))

    def do_g(self):
        # file generated in part F was manually submitted to DAVID
        # and the following enrichment result downloaded
        lines = open('./chart.txt').readlines()[1:]
        terms = list(x.split()[1] for x in lines if len(x.strip()) > 0)
        answer = ("There are {} GO terms enriched in the set of"
                  " differentially regulated genes: {}".format(
                      len(terms), ', '.join(terms)))
        self._write_answer('G', answer)


if __name__ == '__main__':
    ea = ExpressionAnalyzer()
    ea.import_file('../data_set_HL60_U937_NB4_Jurkat.txt')
    for problem in 'abcdefg':
        ea.__getattribute__('do_' + problem)()
