# installing and running Cell Ranger
cd /N/scratch/kmvanden
mkdir tools
cd tools

# https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
# parts of the download command change periodically, copy what is given on the download page
# curl -o *cellranger-7.1.0.tar.gz *https://something*

# unpack the tarball
tar -zxvf cellranger-8.0.1.tar.gz 

# add cell ranger to PATH (only for current session)
cd cellranger-8.0.1
pwd # /N/scratch/kmvanden/tools/cellranger-7.1.0
export PATH="/N/scratch/kmvanden/tools/cellranger-7.1.0:$PATH"
which cellranger # output should be something like /N/scratch/kmvanden/tools/cellranger-7.1.0

# the cellranger count pipeline aligns sequencing reads in FASTQ files to a reference transcriptome and generates a cloupe files
cd mkdir /N/scratch/kmvanden/run_cellranger_count
cd /N/scratch/kmvanden/run_cellranger_count

# download publically available data sets
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
# extract tar file
tar-xvf pbmc_1k_v3_fastqs.tar

# download latest human transcriptome package from 10x website and decompress
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz ​
tar -zxvf refdata-gex-GRCh38-2020-A.tar.gz

# you are now ready to run cellranger count
cellranger count --id=run_count_ikpbmcs --fastqs /N/scratch/kmvanden/run_cellranger_count/pmbc_1k_v3_fastqs/ --sample=pbmc_1k_v3 --transcriptome /N/scratch/kmvanden/run_cellranger_count/refdata-gex-GRCh38-2020-A --create-bam true

# the cell ranger count pipeline outputs are in the outs folder
# within the outs folder there is a filtered_feature_bc_matrix
# matrix.mtx, genes.tsv and barcodes.tsv are in this folder in gzipped format 
# seurat expects Cell Ranger files to be gzipped

# an example script to run cell ranger on slurm
#!/bin/bash​
#SBATCH -A general​
#SBATCH -J pbmc_1k_v3 #name of the job​
#SBATCH -o pbmc_1k_v3_%j.txt #name of the job out-put log file​
#SBATCH -e pbmc_1k_v3_%j.err #name of the job out-put error file​
#SBATCH --mail-type=ALL​
#SBATCH --mail-user=username@iu.edu​
#SBATCH --nodes=1​
#SBATCH --ntasks-per-node=8​
#SBATCH --time=24:00:00​
#SBATCH --partition=general​
#SBATCH --mem=64G​
export PATH=/N/slate/username/tools/cellranger-7.1.0 :$PATH ​
cellranger count --id=run_count_1kpbmcs \ ​
--fastqs /N/slate/username/run_cellranger_count/pbmc_1k_v3_fastqs \​
--sample=pbmc_1k_v3 \ ​
--transcriptome /N/slate/username/run_cellranger_count/refdata-gex-GRCh38-2020-A

