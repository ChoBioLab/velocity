# velocity

Step 1: Setting up a conda env using the requirements file and integrate with jupyter notebook.

conda create --name scvelo --file scvelo_requirements_file.txt
python -m ipykernel install --user --name scvelo --display-name "scvelo_env"

Step 2: Convert Seurat object to anndata

Run SeuratObj_to_anndata.R script to generate python readable files.

Step 3: Run the scvelo python script


