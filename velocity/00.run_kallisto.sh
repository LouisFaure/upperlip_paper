bash scripts/dovelo.sh ~/backup/ML4/fastq/ ML4/velo/ ~/backup/ML4/cellranger/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
bash scripts/dovelo.sh ~/backup/ML5/fastq/ ML5/velo/ ~/backup/ML5/cellranger/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
bash scripts/dovelo.sh ~/backup/ML6/fastq/ ML6/velo/ ~/backup/ML6/cellranger/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
bash scripts/dovelo.sh ~/backup/ML7/fastq/ ML7/velo/ ~/backup/ML7/cellranger/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
bash scripts/dovelo.sh ~/backup/ML8/fastq/ ML8/velo/ ~/backup/ML8/cellranger/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
bash scripts/dovelo.sh ~/backup/ML9/fastq/ ML9/velo/ ~/backup/ML9/cellranger/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
bash scripts/dovelo.sh ~/backup/ML10/fastq/ ML10/velo/ ~/backup/ML10/cellranger/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
bash scripts/dovelo.sh ~/backup/ML11/fastq/ ML11/velo/ ~/backup/ML11/cellranger/outs/filtered_feature_bc_matrix/barcodes.tsv.gz


bash scripts/dovelo_ML4-5.sh ~/backup/ML4/fastq/ ML4/velo/
bash scripts/dovelo_ML4-5.sh ~/backup/ML5/fastq/ ML5/velo/

python scripts/merge_kallisto.py