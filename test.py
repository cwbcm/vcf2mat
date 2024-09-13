import pandas as pd
import numpy as np
from pysam import VariantFile
import re
from vcf import *

ref = pd.read_csv('./refs/ENSEMBL-lite_GRCh38.v94.txt', delimiter='\t', header=0, index_col='gene')



gene_ref = ref.iloc[0]

vcf = VariantFile("../../EpiMu/ukb/VEP.1K.Random.0.reduced.vcf.gz")

samples = vcf.header.samples
dmatrix = pd.DataFrame(np.zeros((len(samples), 1)), index=samples, columns=[gene])
for var in vcf:
    if re.search(r'chr', var.chrom):
        contig = 'chr'+str(gene_ref.chrom)
    else:
        contig = str(gene_ref.chrom)
    break
def _fetch_anno(anno):
    # for fields that could return either direct value or tuple depending on header
    if type(anno) == tuple:
        return anno[0]
    else:
        return anno
for rec in vcf.fetch(contig=contig, start=gene_ref.start, stop=gene_ref.end):
    all_ea = rec.info.get('EA', (None,))
    all_ensp = rec.info.get('Ensembl_proteinid', (rec.info['ENSP'][0],))
    canon_ensp = _fetch_anno(rec.info['ENSP'])
    csq = _fetch_anno(rec.info['Consequence'])
    rec_gene = _fetch_anno(rec.info['SYMBOL'])
    ea = fetch_EA_VEP(all_ea, canon_ensp, all_ensp, csq)
    pass_af_check = af_check(rec, af_field=af_field, max_af=max_af, min_af=min_af)
    if not np.isnan(ea).all() and gene == rec_gene and pass_af_check:
        gts = pd.Series([convert_zygo(rec.samples[sample]['GT']) for sample in samples], index=samples, dtype=int)
        dmatrix[gene] += ea*gts

    def _fetch_anno(anno):
        # for fields that could return either direct value or tuple depending on header
        if type(anno) == tuple:
            return anno[0]
        else:
            return anno

def fetch_EA_VEP(EA, canon_ensp, all_ensp, csq):
    if 'stop_gained' in csq or 'frameshift_variant' in csq or 'stop_lost' in csq:# or 'splice_donor_variant' in csq or 'splice_acceptor_variant' in csq:
        return 100
    try:
        canon_idx = all_ensp.index(canon_ensp)
    except ValueError:
        return np.nan
    else:
        return validate_EA(EA[canon_idx])

def af_check(rec, af_field, max_af, min_af):
    """
    Check if variant allele frequency passes filters
    Args:
        rec (VariantRecord)
        af_field (str): Name of INFO field containing allele frequency information
        max_af (float): Maximum allele frequency for variant
        min_af (float): Minimum allele frequency for variant
    Returns:
        bool: True of AF passes filters, otherwise False
    """

    af = rec.info['AF']
    if type(af) == tuple:
        af = af[0]
    try:
        return min_af < float(af) < max_af
    except:
        return False

def validate_EA(ea):
    """
    Checks for valid EA score
    Args:
        ea (str/float/None): EA score as string
    Returns:
        float: EA score between 0-100 if valid, otherwise returns NaN
    """
    try:
        ea = float(ea)
    except ValueError:
        if type(ea) == str and (ea == 'fs-indel' or 'STOP' in ea):
            ea = 100
        else:
            ea = np.nan
    except TypeError:
        ea = np.nan
    return ea

def convert_zygo(genotype):
    """
    Convert a genotype tuple to a zygosity integer
    Args:
        genotype (tuple): The genotype of a variant for a sample
    Returns:
        int: The zygosity of the variant (0/1/2)
    """
    if genotype in [(1, 0), (0, 1)]:
        zygo = 1
    elif genotype == (1, 1):
        zygo = 2
    else:
        zygo = 0
    return zygo

gene = "OR2T10"
sumEA_mat = pd.DataFrame(np.zeros((len(samples), 1)), index=samples, columns=[gene])
maxEA_mat = pd.DataFrame(np.zeros((len(samples), 1)), index=samples, columns=[gene])
for rec in vcf.fetch(contig=contig, start=gene_ref.start, stop=gene_ref.end):
    all_ea = rec.info.get('EA', (None,))
    all_ensp = rec.info.get('Ensembl_proteinid', (rec.info['ENSP'][0],))
    canon_ensp = _fetch_anno(rec.info['ENSP'])
    csq = _fetch_anno(rec.info['Consequence'])
    rec_gene = _fetch_anno(rec.info['SYMBOL'])
    ea = fetch_EA_VEP(all_ea, canon_ensp, all_ensp, csq)
    pass_af_check = af_check(rec, af_field="AF", max_af=1, min_af=0)
    if not np.isnan(ea).all() and pass_af_check:
        gts = pd.Series([convert_zygo(rec.samples[sample]['GT']) for sample in samples], index=samples, dtype=int)
        maxEA_mat[gene] = np.maximum(maxEA_mat[gene], ea*(gts>0))
        sumEA_mat[gene] += ea*gts
    print(rec)


mat1, mat2 = vep_to_mat("../../EpiMu/ukb/VEP.1K.Random.0.reduced.vcf.gz",
                        gene = "OR2T10", gene_ref = ref.iloc[0],
                        samples = samples, max_af = 1, min_af = 0, af_field='AF')
