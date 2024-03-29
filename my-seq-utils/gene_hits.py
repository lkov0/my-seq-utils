#read bam file and list of regions, output regions with hits, percentage of regions with hits
# usage: python gene_hits.py bamfile bedfile

from sys import argv
import pysam
import pandas as pd
import numpy as np

script, bam, bed = argv

def main():
    bedstem = bed.replace('.bed', '')
    bamstem = bam.replace('.bam', '')
    bamfile = pysam.AlignmentFile(bam, 'rb')
    regions = pd.read_csv(bed, delimiter = "\t", header = None)

    regions['coverage'] = np.vectorize(getcov)(bamfile, regions[0], regions[1], regions[2])

    print(f"Writing bed file with coverage to: {bamstem + '.coverage.bed'}")
    writecov(bamstem + 'coverage.bed', regions)

    print(f"Proportion of regions with hits: {covperc(regions)}")


# returns number of reads that aligned to the specified region
def getcov(df, chr, sta, sto):
    return df.count(contig = chr, start = sta, stop = sto)

# calculates the proportion of regions with hits
def covperc(df):
    expr_regions = len(df[(df['coverage'] > 0)])
    return expr_regions / len(df)

# writes bedfile with coverage information to a tsv file
def writecov(fname, df):
    df = df[[0, 1, 2, 'coverage']]
    df.columns = ['chrom', 'start', 'stop', 'coverage']
    df.to_csv(fname, sep = "\t", index = False)

if __name__=="__main__":
   main()
