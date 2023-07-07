# velocity

## Step 0: Generate loom files:
- Run the following on terminal:
  
ml samtools/1.17

ml velocyto

source activate /hpc/packages/minerva-centos7/velocyto/0.17/velocyto

### can either run directly on terminal or as an interactive job
bsub -P acc_untreatedIBD -q interactive -n 6 -R "span[hosts=1]" -R rusage[mem=8000] -W 3:00 -Is /bin/bash
### run the code below - this will generate loom files in the same directry as the cellranger samples are located
velocyto run10x </path/to/cellranger/sample/directory> <path/to/the/reference/genome>/genes.gtf

# example:

velocyto run10x /sc/arion/projects/untreatedIBD/stalware/Josh_csf2rb/analysis/cellranger/cr_count_2022-05-31_1423/J2/outs/filtered_feature_bc_matrix /sc/arion/projects/untreatedIBD/stalware/ref_genome/Danio.rerio_genome/genes.gtf


## Step 1: Setting up a conda env using the requirements file and integrate with jupyter notebook.

conda create --name scvelo --file scvelo_requirements_file.txt

python -m ipykernel install --user --name scvelo --display-name "scvelo_env"

## Step 2: Convert Seurat object to anndata

Run SeuratObj_to_anndata.R script to generate python readable files.

Step 3: Run the scvelo script

