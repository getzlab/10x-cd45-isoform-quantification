## Download data from https://www.10xgenomics.com/resources/datasets/pbm-cs-of-a-healthy-donor-5-gene-expression-and-cell-surface-protein-1-standard-3-0-0
## PBMCs of a Healthy Donor - 5' Gene Expression with a Panel of TotalSeq™-C Antibodies

CD45er_path="..." ## Add you CD45er path

wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-vdj/3.0.0/vdj_v1_hs_pbmc2_5gex_protein/vdj_v1_hs_pbmc2_5gex_protein_possorted_genome_bam.bam
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-vdj/3.0.0/vdj_v1_hs_pbmc2_5gex_protein/vdj_v1_hs_pbmc2_5gex_protein_possorted_genome_bam.bam.bai
wget https://cf.10xgenomics.com/samples/cell-vdj/3.0.0/vdj_v1_hs_pbmc2_5gex_protein/vdj_v1_hs_pbmc2_5gex_protein_filtered_feature_bc_matrix.h5

python $CD45er_path/run.py \
    vdj_v1_hs_pbmc2_5gex_protein_possorted_genome_bam.bam \
    $CD45er_path/reference/CD45_exons_nochr.hg38.txt

