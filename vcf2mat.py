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
    parser.add_argument('--ref', nargs='?', default='hg38', help='genome reference file. Use hg38 or hg19 for saved ref, otherwise, provide path.')
    parser.add_argument('--savepath', nargs='?', default='./',help='save path for output')
    parser.add_argument('--minaf', type=float, default=0, help='minimum allele frequency cutoff')
    parser.add_argument('--maxaf', type=float, help='maximum allele frequency cutoff')
    parser.add_argument('--Ann',default='VEP', choices=('VEP', 'ANNOVAR'),help='EA annotation method')
    parser.add_argument('--cores', type=int, default=1, help='number of CPUs to use for multiprocessing')
    parser.add_argument('--prefix', nargs='?', default='',help='prefix to output files')
    parser.add_argument('--output', nargs='?', default='pEA,sumEA,maxEA,silent', help='matrices to output, separate by ,')
    parser.add_argument('--format', default='tsv', choices=('tsv', 'parquet'), help='output format of the matrices')
    return parser.parse_args()



def main(args):
    if args.Ann=='ANNOVAR':
        if args.ref=='hg19':
            ref = pd.read_csv('./refs/refGene-lite_hg19.May2013.txt', delimiter='\t', header=0, index_col='gene')
        elif args.ref=='hg38':
            ref = pd.read_csv('./refs/refGene-lite_hg38.June2017.txt', delimiter='\t', header=0, index_col='gene')
        else:
            ref = pd.read_csv(args.ref, delimiter='\t', header=0, index_col='gene')
    elif args.Ann=='VEP':
        if args.ref=='hg19':
            ref = pd.read_csv('./refs/ENSEMBL-lite_GRCh37.v75.txt', delimiter='\t', header=0, index_col='gene')
        elif args.ref=='hg38':
            ref = pd.read_csv('./refs/ENSEMBL-lite_GRCh38.v94.txt', delimiter='\t', header=0, index_col='gene')
        else:
            ref = pd.read_csv(args.ref, delimiter='\t', header=0, index_col='gene')
    
    # samples = pd.read_csv(args.samples, header=None, index_col=0)
    # total_samples = samples.index.astype(str).tolist()
    total_samples = list(VariantFile(args.VCF).header.samples)
    # if args.Ann=='ANNOVAR':
    #    results = Parallel(n_jobs=args.cores)(delayed(parse_ANNOVAR)(args.VCF, gene, ref.loc[gene], total_samples, min_af=0, max_af=args.maxaf,af_field='AF',EA_parser='canonical') for gene in tqdm(ref.index.unique()))
    if args.Ann=='VEP':
        results = Parallel(n_jobs=args.cores)(delayed(vep_to_mat)(args.VCF, gene, ref.loc[gene], total_samples, max_af=args.maxaf, min_af=args.minaf) for gene in tqdm(ref.index.unique()))

    sumEA_mat = pd.concat([result[0] for result in results], axis=1)
    maxEA_mat = pd.concat([result[1] for result in results], axis=1)
    pEA_mat = pd.concat([result[2] for result in results], axis=1)
    silent_mat = pd.concat([result[3] for result in results], axis=1)
    
    if 'pEA' in args.output.split(","):
        pEA_pl = pl.from_pandas(pEA_mat)
        pEA_pl = pEA_pl.with_columns(pl.all().round(4))
        pEA_pl.insert_column(0, pl.from_pandas(pEA_mat.index).alias("sample"))
        if args.format == "tsv":
            pEA_pl.write_csv(args.savepath+f'/{args.prefix}pEA_mat.tsv', separator='\t', float_precision=4)
        else:
            pEA_pl.write_parquet(args.savepath+f'/{args.prefix}pEA_mat.parquet')
            
    if 'sumEA' in args.output.split(","):
        sumEA_pl = pl.from_pandas(sumEA_mat)
        sumEA_pl = sumEA_pl.with_columns(pl.all().round(4))
        sumEA_pl.insert_column(0, pl.from_pandas(sumEA_mat.index).alias("sample"))
        if args.format == "tsv":
            sumEA_pl.write_csv(args.savepath+f'/{args.prefix}sumEA_mat.tsv', separator='\t', float_precision=4)
        else:
            sumEA_pl.write_parquet(args.savepath+f'/{args.prefix}sumEA_mat.parquet')

    if 'maxEA' in args.output.split(","):
        maxEA_pl = pl.from_pandas(maxEA_mat)
        maxEA_pl = maxEA_pl.with_columns(pl.all().round(4))
        maxEA_pl.insert_column(0, pl.from_pandas(maxEA_mat.index).alias("sample"))
        if args.format == "tsv":
            maxEA_pl.write_csv(args.savepath+f'/{args.prefix}maxEA_mat.tsv', separator='\t', float_precision=4)
        else:
            maxEA_pl.write_parquet(args.savepath+f'/{args.prefix}maxEA_mat.parquet')

    if 'silent' in args.output.split(","):
        silent_pl = pl.from_pandas(silent_mat)
        silent_pl = silent_pl.with_columns(pl.all().round(4))
        silent_pl.insert_column(0, pl.from_pandas(silent_mat.index).alias("sample"))
        if args.format == "tsv":
            silent_pl.write_csv(args.savepath+f'/{args.prefix}silent_mat.tsv', separator='\t', float_precision=4)
        else:
            silent_pl.write_parquet(args.savepath+f'/{args.prefix}silent_mat.parquet')

if __name__ == "__main__":
    args = parse_args()
    main(args)
