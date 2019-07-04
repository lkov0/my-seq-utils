#!/usr/bin/env python3

# This script takes a DAGChainer output file from SynMap (CoGe) as input and outputs a list of syntenic gene pairs.
# DAGChainer table should have columns 'chr1 gene1 start1 stop1 chr2 gene2 start2 stop2 eValue score'
# genePair table should be a tsv file with genes separated by '\t', one pair per line
# 'start' and 'stop' values should be reflective of gene order and not genomic coordinates

import argparse
import pandas as pd
import bisect
from collections import defaultdict


class sortable_dag_row (pd.Series):
    """ class for DAG rows which can be sorted on start2.
        supports comparison by either integer or row; in latter case, start2 is extracted for comparison
    """

    def __init__(self, s):
        super().__init__(s)

    def __lt__(self, o):
        if not isinstance(o, int):
            o = int(o['start2'])
        return int(self['start2']) < o

    def __gt__(self, o):
        if not isinstance(o, int):
            o = int(o['start2'])
        return int(self['start2']) > o

    def __eq__(self, o):
        if not isinstance(o, int):
            o = int(o['start2'])
        return int(self['start2']) == o


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--synfile", required=True,
                        help="DAGChainer output from CoGe SynMap")
    parser.add_argument("-p", "--pairfile", required=True,
                        help="tsv file with one gene pair per line, genes separated by '\t'")
    parser.add_argument("-o", "--outputDir", required=True,
                        help="Desired output directory")
    args = parser.parse_args()

    dagnames = ['chr1', 'gene1', 'start1', 'stop1', 'chr2',
                'gene2', 'start2', 'stop2', 'eValue', 'score']

    dag = pd.read_csv(args.synfile, sep='\t', comment='#', names=dagnames)
    genepairs = pd.read_csv(args.pairfile, sep='\t', names=['g1', 'g2'])

    # create an index that is essentially { gene: { chr: [...] } }, where the inner lists are sorted on start2
    binned = defaultdict(lambda: defaultdict(lambda: []))
    for _, d in dag.iterrows():
        bisect.insort_left(binned[str(d['gene1'])]
                           [str(d['chr2'])], sortable_dag_row(d))

    syn_pairs = pd.DataFrame(
        columns=['g1', 'g2', 's1', 's2'], data=find_all_syntenic_pairs(genepairs, binned))
    syn_pairs.to_csv(args.outputDir, sep="\t", index=False)


def find_all_syntenic_pairs(genepairs, binned):
    """ generates rows containing syntenic pairs
    """
    for _, row in genepairs.iterrows():
        yield from find_syntenic_pairs(binned, str(row['g1']), str(row['g2']))


def find_syntenic_pairs(binned, g1, g2):
    """ generates rows containing gene pairs syntenic to specified pair
    """

    g1_syn = binned[g1]
    g2_syn = binned[g2]

    for chrom, rows_g1 in g1_syn.items():
        rows_g2 = g2_syn[chrom]
        for row_g1 in rows_g1:

            # remember the index where we found (start2 - 1) so we don't have to search for (start2 + 1) from the beginning
            i = 0

            for start2 in [int(row_g1['start2']) - 1, int(row_g1['start2']) + 1]:

                # perform binary search for desired start2 value
                i = bisect.bisect_left(rows_g2, start2)

                # when there is no match, bisect may return an index either (a) beyond the list bounds, or (b) preceding our desired value
                if i < len(rows_g2) and rows_g2[i] == start2:
                    row_g2 = rows_g2[i]
                    yield [row_g1['gene1'], row_g2['gene1'], row_g1['gene2'], row_g2['gene2']]


if __name__ == '__main__':
    main()
