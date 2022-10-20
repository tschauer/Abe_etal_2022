
### copy these to the shell ###


# conda install mamba -n base -c conda-forge
#
# mamba env create --name ken_rnaseq --file ken_rnaseq.yaml
# mamba env update --name ken_rnaseq --file ken_rnaseq.yaml


conda activate ken_rnaseq

### R packages from github ###
# R
# library(devtools)
# install_github("tschauer/HelpersforDESeq2")
# install_github("js229/Vennerable")


################################################################################


snakemake -np > snakemake_log_full.txt
#snakemake --rulegraph | dot -Tpng > snakemake_rulegraph.png


sbatch --partition=normal_q \
--time=24:00:00  \
--constraint=Lustre_File_System \
--wrap="snakemake -j 30 --cluster-config cluster.json --cluster 'sbatch --constraint=Lustre_File_System -p {cluster.partition} --exclude {cluster.exclude} -n {cluster.n} --mem {cluster.mem} -t {cluster.time}'"


################################################################################
