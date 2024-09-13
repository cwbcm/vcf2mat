import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm
import argparse
from vcf import *
import polars as pl


def parse_args():
    """
    Parses the arguments.
    """
    parser = argparse.ArgumentParser(description="Mu arguments")
    parser.add_argument('--VCF', nargs='?', default='./',help='Location of Cohort VCF')
    parser.add_argument('--samples', nargs='?', default='./samples.txt',help='samples file path')
    parser.add_argument('--GeneLength', nargs='?', default='./refs/gene_length.csv',help='gene length file path')
    parser.add_argument('--ref', nargs='?', default='hg38', choices=('hg19', 'hg38'), help='genome reference file')
    parser.add_argument('--savepath', nargs='?', default='./',help='save path for output')
    parser.add_argument('--minaf', type=float, default=0, help='minimum allele frequency cutoff')
    parser.add_argument('--maxaf', type=float, help='maximum allele frequency cutoff')
    parser.add_argument('--Ann',default='VEP', choices=('VEP', 'ANNOVAR'),help='EA annotation method')
    parser.add_argument('--cores', type=int, default=1, help='number of CPUs to use for multiprocessing')
    parser.add_argument('--chrX', type=int,default=1, help='1 if there is sex chromosome in the VCF, 0 if there is no sex chromosome in the VCF')
    parser.add_argument('--prefix', nargs='?', default='',help='prefix to output files')
    parser.add_argument('--output', nargs='?', default='pEA,sumEA,maxEA', help='matrices to output, separate by ,')
    return parser.parse_args()



def main(args):
    if args.chrX==1:
        if args.Ann=='ANNOVAR':
            if args.ref=='hg19':
                ref = pd.read_csv('./refs/refGene-lite_hg19.May2013.txt', delimiter='\t', header=0, index_col='gene')
            elif args.ref=='hg38':
                ref = pd.read_csv('./refs/refGene-lite_hg38.June2017.txt', delimiter='\t', header=0, index_col='gene')
        elif args.Ann=='VEP':
            if args.ref=='hg19':
                ref = pd.read_csv('./refs/ENSEMBL-lite_GRCh37.v75.txt', delimiter='\t', header=0, index_col='gene')
            elif args.ref=='hg38':
                ref = pd.read_csv('./refs/ENSEMBL-lite_GRCh38.v94.txt', delimiter='\t', header=0, index_col='gene')
    elif (args.chrX==0) & (args.ref=='hg38') & (args.Ann=='VEP'):
        ref = pd.read_csv('./refs/ENSEMBL-lite_GRCh38.v94.noX.txt', delimiter='\t', header=0, index_col='gene')

    samples = pd.read_csv(args.samples, header=None, index_col=0)
    total_samples = samples.index.astype(str).tolist()
    # if args.Ann=='ANNOVAR':
    #    results = Parallel(n_jobs=args.cores)(delayed(parse_ANNOVAR)(args.VCF, gene, ref.loc[gene], total_samples, min_af=0, max_af=args.maxaf,af_field='AF',EA_parser='canonical') for gene in tqdm(ref.index.unique()))
    if args.Ann=='VEP':
        results = Parallel(n_jobs=args.cores)(delayed(vep_to_mat)(args.VCF, gene, ref.loc[gene], total_samples, max_af=args.maxaf, min_af=args.minaf) for gene in tqdm(ref.index.unique()))

    sumEA_mat = pd.concat([result[0] for result in results], axis=1)
    maxEA_mat = pd.concat([result[1] for result in results], axis=1)
    pEA_mat = pd.concat([result[2] for result in results], axis=1)
    
    if 'pEA' in args.output.split(","):
        pEA_pl = pl.from_pandas(pEA_mat)
        pEA_pl.insert_column(0, pl.from_pandas(pEA_mat.index).alias("sample"))
        pEA_pl.write_csv(args.savepath+f'/{args.prefix}pEA_mat.tsv', separator='\t', float_precision=4)
    if 'sumEA' in args.output.split(","):
        sumEA_pl = pl.from_pandas(sumEA_mat)
        sumEA_pl.insert_column(0, pl.from_pandas(sumEA_mat.index).alias("sample"))
        sumEA_pl.write_csv(args.savepath+f'/{args.prefix}sumEA_mat.tsv', separator='\t', float_precision=2)
    if 'maxEA' in args.output.split(","):
        maxEA_pl = pl.from_pandas(maxEA_mat)
        maxEA_pl.insert_column(0, pl.from_pandas(maxEA_mat.index).alias("sample"))
        maxEA_pl.write_csv(args.savepath+f'/{args.prefix}maxEA_mat.tsv', separator='\t', float_precision=2)

if __name__ == "__main__":
    args = parse_args()
    main(args)
