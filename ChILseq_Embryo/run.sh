### copy these to the shell ###

# mamba env create --name ken_chilseq --file ../env_ken_chilseq.yaml
# mamba env update --name ken_chilseq --file ../env_ken_chilseq.yaml


conda activate ken_chilseq

# conda deactivate
# conda env remove --name ken_chilseq



################################################################################


snakemake -np > snakemake_log.txt
#snakemake --rulegraph | dot -Tpng > snakemake_rulegraph.png


sbatch --partition=normal_q \
--time=24:00:00  \
--constraint=Lustre_File_System \
--wrap="snakemake -j 50 --cluster-config cluster.json --cluster 'sbatch --constraint=Lustre_File_System -p {cluster.partition} -q {cluster.q} --exclude {cluster.exclude} -n {cluster.n} --mem {cluster.mem} -t {cluster.time}'"


################################################################################
