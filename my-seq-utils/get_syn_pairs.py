# This script takes a DAGChainer output file from SynMap (CoGe) as input and outputs a list of syntenic gene pairs.
# DAGChainer table should have columns 'chr1 gene1 start1 stop1 chr2 gene2 start2 stop2 eValue score'
# genePair table should be a tsv file with genes separated by '\t', one pair per line
# 'start' and 'stop' values should be reflective of gene order and not genomic coordinates

import argparse
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument("-s", "--synfile", help = "DAGChainer output from CoGe SynMap")
parser.add_argument("-p", "--pairfile", help = "tsv file with one gene pair per line, genes separated by '\t'")
parser.add_argument("-o", "--outputDir", help = "Desired output directory")
args = parser.parse_args()


def main():

    dagnames = ['chr1','gene1','start1','stop1','chr2','gene2','start2','stop2','eValue','score']

    dag = pd.read_csv(args.synfile, sep='\t', comment='#', names=dagnames)
    genepairs = pd.read_csv(args.pairfile, sep='\t', names=['g1', 'g2'])

    syn_pairs = pd.DataFrame(columns=['g1', 'g2', 's1', 's2'])

    count = 0;
    for i, row in genepairs.iterrows():
        count = count + 1
        print(f"ran find function {count} times")
        syn_temp = find_syntenic_pairs(dag, row['g1'], row['g2'])
        syn_pairs = pd.concat([syn_pairs, syn_temp])

    syn_pairs.to_csv(args.output, sep="\t", index=False)

# returns a dataframe containing gene pairs syntenic to specified pairs


def find_syntenic_pairs(df, g1, g2):

    # fetch syntenic pairs for gene1 and gene2
    g1_syn = df.loc[df['gene1'] == str(g1)]
    g2_syn = df.loc[df['gene1'] == str(g2)]

    df_syn = pd.DataFrame(columns=['g1', 'g2', 's1', 's2'])

    if not g1_syn.empty:
        print("running synteny find")
        for x, row_g1 in g1_syn.iterrows():
            for y, row_g2 in g2_syn.iterrows():
                # check if gene1 and gene2 contain syntenic pairs on same chromosome
                if is_syntenic_pair(row_g1, row_g2):
                    print("syntenic pairs found!")
                    new_row = [row_g1['gene1'], row_g2['gene1'], row_g1['gene2'], row_g2['gene2']]
                    df_syn.loc[len(df_syn)] = new_row

    return df_syn

# checks if a syntenic match exists between datapoint inputs


def is_syntenic_pair(row_g1, row_g2):
    return (row_g1['chr2'] == row_g2['chr2']  # are syntenic matches on the same chromosome?
            and (int(row_g1['start2']) == int(row_g2['start2']) + 1  # are syntenic matches adjacent on that chromosome?
                or int(row_g1['start2']) == int(row_g2['start2']) - 1))


if __name__ == '__main__':
    main()
