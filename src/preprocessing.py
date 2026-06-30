import allel
import pandas as pd
import numpy as np

vkorc1 = allel.read_vcf('/mnt/c/Users/user/genomic-variant-ml/data/raw/vkorc1_rs9923231.vcf.gz')
cyp2c9 = allel.read_vcf('/mnt/c/Users/user/genomic-variant-ml/data/raw/cyp2c9_variants.vcf.gz')

vkorc1_gt = allel.GenotypeArray(vkorc1['calldata/GT']) #the genotype field for every sample at every variant
cyp2c9_gt = allel.GenotypeArray(cyp2c9['calldata/GT'])

vkorc1_samples = vkorc1['samples']
cyp2c9_samples = cyp2c9['samples']

#counts ALT alleles per sample per cariant: 0, 1, or 2
vkorc1_dosage = vkorc1_gt.to_n_alt()

#building vkorc1 genotype lookup
def get_vkorc1_genotype(gt_call, ref='G', alt='A'):
    """
    Convert a 0/0, 0/1, 1/1 call into a GG/GA/AA string.
    """
    alleles = sorted([ref if a == 0 else alt for a in gt_call], reverse=True)
    return ''.join(alleles)

vkorc1_results = []
for i, sample_id in enumerate(vkorc1_samples):
    genotype_str = get_vkorc1_genotype(vkorc1_gt[0, i])
    vkorc1_results.append({'sample_id': sample_id, 'vkorc1_genotype': genotype_str})

vkorc1_df = pd.DataFrame(vkorc1_results)

#building cyp2c9 stae-allele table
def get_cyp2c9_star_allele(rs1799853_gt, rs1057910_gt):
    """
    Reconstruct cyp2c9 star-allele genotype (e.g. '*1/*1', '*1/*2', '*2/*3')
    Parameters
    ----------
    rs1799853_gt : tuple/list of two ints, e.g. (0, 0), (0, 1), (1, 1)
        Genotype call at rs1799853 (defines *2 allele).
        0 = reference allele, 1 = alt allele (*2)
    rs1057910_gt : tuple/list of two ints, e.g. (0, 0), (0, 1), (1, 1)
        Genotype call at rs1057910 (defines *3 allele).
        0 = reference allele, 1 = alt allele (*3)

    Returns
    -------
    str
        Star-allele genotype string, e.g. '*1/*1', '*1/*2', '*2/*3'
        Returns None if the combination is biologically inconsistent
        (e.g. both *2 and *3 alt alleles present in a way that exceeds 2 total alleles)
    """

    n_star2 = sum(rs1799853_gt)
    n_star3 = sum(rs1057910_gt)

    total_variant_alleles = n_star2 + n_star3

    if total_variant_alleles > 2:
        return None
    
    # num of normal function (*1) alleles = whatever's left over
    n_star1 = 2 - total_variant_alleles

    alleles = (['*1'] * n_star1) + (['*2'] * n_star2) + (['*3'] * n_star3)
    alleles_sorted = sorted(alleles, key=lambda a: int(a[1]))

    return f"{alleles_sorted[0]}/{alleles_sorted[1]}"

cyp2c9_results = []
for i, sample_id in enumerate(cyp2c9_samples):
    rs1799853_gt = cyp2c9_gt[0, i]
    rs1057910_gt = cyp2c9_gt[1, i]
    star_allele = get_cyp2c9_star_allele(rs1799853_gt, rs1057910_gt)
    cyp2c9_results.append({'sample_id': sample_id, 'cyp2c9_genotype': star_allele})

cyp2c9_df = pd.DataFrame(cyp2c9_results)

#merging vkorc1 + cyp2c9 + dosing table to generate labels
combined_df = vkorc1_df.merge(cyp2c9_df, on='sample_id')

dosing_table = pd.read_csv('/mnt/c/Users/user/genomic-variant-ml/data/raw/warfarin_dosing_table.csv')

#assign sensitivity_class based on genotype combination
combined_df = combined_df.merge(dosing_table, left_on=['vkorc1_genotype', 'cyp2c9_genotype'],
                                right_on=['vkorc1_genotype', 'cyp2c9_genotype'], how='left')

print("Unmatched rows:", combined_df['sensitivity_class'].isna().sum())
combined_df = combined_df.dropna(subset=['sensitivity_class'])

feature_rows = []
for i, sample_id in enumerate(vkorc1_samples):
    feature_rows.append({'sample_id': sample_id, 'vkorc1_dosage': sum(vkorc1_gt[0, i]),
                         'cyp2c9_2_dosage': sum(cyp2c9_gt[0, i]), 'cyp2c9_3_dosage': sum(cyp2c9_gt[1, i])})

feature_df = pd.DataFrame(feature_rows)
final_df = feature_df.merge(combined_df[['sample_id', 'sensitivity_class']], on='sample_id')
final_df.to_csv('/mnt/c/Users/user/genomic-variant-ml/data/processed/feature_matrix.csv', index=False)

from sklearn.model_selection import train_test_split

X = final_df[['vkorc1_dosage', 'cyp2c9_2_dosage', 'cyp2c9_3_dosage']]
y = final_df['sensitivity_class']

X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, test_size=0.2, random_state=42)

print(final_df['sensitivity_class'].value_counts())