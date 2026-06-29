<!-- cd data/raw

# VKORC1 — chromosome 16
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi

# CYP2C9 — chromosome 10 (confirm version suffix in the directory listing first)
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi

sudo apt install bcftools   # if not already installed

bcftools view -r 16:31100000-31110000 ALL.chr16.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz -Oz -o vkorc1_region.vcf.gz

bcftools view -r 10:96695000-96750000 ALL.chr10.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz -Oz -o cyp2c9_region.vcf.gz

bcftools view -i 'ID="rs9923231"' vkorc1_region.vcf.gz -Oz -o vkorc1_rs9923231.vcf.gz
bcftools view -i 'ID="rs1799853" || ID="rs1057910"' cyp2c9_region.vcf.gz -Oz -o cyp2c9_variants.vcf.gz

wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel -->