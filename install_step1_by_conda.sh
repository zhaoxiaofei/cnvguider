#!/usr/bin/env bash

if [ -z "${1}" ] ; then
    envname=cnvguider;
else 
    envname=$1;
fi
envparam="-n ${envname}"

bioconda_packages="hmmcopy bioconductor-hmmcopy bioconductor-ctc bioconductor-dnacopy bioconductor-copynumber bioconductor-wgsmapp bioconductor-scope bioconductor-helloranges"
condaforge_packages="apache-airflow docopt picard pyfaidx r-devtools r-inline r-gplots r-scales r-plyr r-ggplot2 r-gridExtra r-fastcluster r-heatmap3"

conda=mamba
condaforge="" #" -c conda-forge " # conda-forge
bioconda="" #"-c conda-forge -c bioconda"

$conda create -n $envname "python<3.13"  matplotlib numpy 'pandas<2.0' scipy scikit-learn seaborn tasklogger vim bcftools bedtools bowtie2 bwa gatk4 htslib samtools htslib pysam ${condaforge_packages} ${bioconda_packages}

#$conda install --yes $envparam $condaforge 
#$conda install --yes $envparam $bioconda bcftools bedtools bowtie2 bwa gatk4 htslib samtools htslib # pysam 

# HMMcopy CopyNumber
#$conda install --yes $envparam $bioconda bioconductor-hmmcopy hmmcopy bioconductor-copynumber

# Ginkgo
#$conda install --yes $envparam $condaforge r-devtools

# https://github.com/deepomicslab/SeCNV
#$conda install --yes $envparam docopt picard pyfaidx

#$conda install --yes $envparam \
#    bioconductor-ctc bioconductor-dnacopy bioconductor-wgsmapp bioconductor-scope bioconductor-helloranges \
#    r-inline r-gplots r-scales r-plyr r-ggplot2 r-gridExtra r-fastcluster r-heatmap3

# https://github.com/biosinodx/SCCNV 
# echo nothing new to be installed for SCCNV

# SCOPE
# echo nothing new to be installed for SCOPE # $conda install --yes $envparam $bioconda bioconductor-scope # install with R to avoid error

# https://github.com/xikanfeng2/SCYN
$conda install -n $envname --yes tasklogger
$conda run     -n $envname       pip install scyn

#sed -i 's;all_cnv = all_cnv.append(cnv);all_cnv = pd.concat([all_cnv, cnv]);g' ${CONDA_PREFIX}/lib/python3.*/site-packages/scyn/utils.py
#sed -i 's;all_cnv = all_cnv.append(cnv);all_cnv = pd.concat([all_cnv, cnv]);g' ${CONDA_PREFIX}/lib/python3.*/site-packages/scyn/utils.py

if false; then
    conda list       -n ${envname} -e | grep -v "^scyn="   > env/requirements.list_e_no_pypi.txt
    conda env export -n ${envname}                         > env/freeze.env_export.yml
    conda env export -n ${envname} --no-builds             > env/freeze.env_export_no_builds.yml
    conda env export -n ${envname} --from-history          > env/freeze.env_export_from_history.yml
fi

conda activate $envname && sed -i 's;Gini<=0.12;Gini<=0.21;g' ${CONDA_PREFIX}/lib/python3.*/site-packages/scyn/utils.py
conda activate $envname && sed -i 's;perform_qc(Y_raw = Y_raw,;perform_qc(mapq20_thresh = 0.1, Y_raw = Y_raw,;g' ${CONDA_PREFIX}/lib/python3.*/site-packages/scyn/utils.py

exit 0



script_dir=$(dirname "$(realpath "$0")")
#/home/zhaoxiaofei/miniconda3/envs/test/lib/R/library/HMMcopy/doc/HMMcopy.R

# https://github.com/shahcompbio/single_cell_pipeline

pushd "${script_dir}/../software"
wget https://github.com/shahcompbio/single_cell_pipeline/archive/refs/tags/v0.8.26.tar.gz
tar -xvf v0.8.26.tar.gz
popd

# conda install -c shahcompbio pypeliner

#pushd /mnt/d/software/
#rm 1.0.tar.gz || true
wget https://github.com/KChen-lab/MEDALT/archive/refs/tags/1.0.tar.gz
tar -xvf 1.0.tar.gz
#git clone https://github.com/KChen-lab/MEDALT.git # 3d4a6d548171ede333310d2ef25c12cdccd11a2b
Rscript -e 'install.packages("igraph")'
#Rscript -e 'BiocManager::install("HelloRanges")'



# CONDA_PREFIX==/home/zhaoxiaofei/miniconda3/envs/single_cell_pipeline-0.8.26/

# tree softwares

wget https://github.com/shahcompbio/single_cell_pipeline/archive/refs/tags/v0.8.26.tar.gz
tar -xvf v0.8.26.tar.gz

#rm 1.0.tar.gz || true
#wget https://github.com/KChen-lab/MEDALT/archive/refs/tags/1.0.tar.gz
#tar -xvf 1.0.tar.gz
#git clone https://github.com/KChen-lab/MEDALT.git # 3d4a6d548171ede333310d2ef25c12cdccd11a2b
#Rscript -e 'install.packages("igraph")'
#Rscript -e 'BiocManager::install("HelloRanges")'

popd

