# Ginkgo

#git clone https://github.com/robertaboukhalil/ginkgo.git && pushd ginkgo && git checkout d7c7790 && make && popd

sed 's;library(gridExtra);library(gridExtra)\n\n#NOTE: The workaround https://github.com/ChristophH/gplots described at https://github.com/robertaboukhalil/ginkgo no longer works with R version 4.X.X\n#NOTE: Therefore, we fall back to the heatmap function\nheatmap.2 = function(...) { return(heatmap(...)) };g' ginkgo/scripts/process.R -i

sed 's;tar -czf;#NOTE: the tar command results in runtime error and is therefore commented out\n#tar -czf;g' ginkgo/cli/ginkgo.sh -i

mkdir -p ../../refs/refs3to4/
pushd    ../../refs/refs3to4/
#wget https://labshare.cshl.edu/shares/schatzlab/www-data/ginkgo/genomes/hg19.tgz
#tar -xvf hg19.tgz
popd
rm -r ${PWD}/ginkgo/genomes/hg19 || true && cp -rs ${PWD}/../../refs/refs3to4/hg19 ${PWD}/ginkgo/genomes/ || true

# HMMCopy

wget https://github.com/shahcompbio/single_cell_pipeline/archive/refs/tags/v0.8.26.tar.gz
tar -xvf v0.8.26.tar.gz # use sudo tar in case of permission error
conda create --yes --name single_cell_pipeline-0.8.26 --file single_cell_pipeline-0.8.26/docker/hmmcopy/conda_requirements.yml
conda create --yes --name hmmcopy --file hmmcopy_919ca93_requirements.yml # from https://github.com/mondrian-scwgs/mondrian/blob/main/docker/hmmcopy/requirements.yml
git clone https://github.com/shahcompbio/hmmcopy_utils.git
mkdir -p hmmcopy_utils/build && pushd hmmcopy_utils/build && cmake .. && make -j 8 && popd && cp -r hmmcopy_utils/build/bin hmmcopy_utils/ # use sudo cmake if make generates permission error 

# https://github.com/biosinodx/SCCNV
git clone https://github.com/biosinodx/SCCNV.git && pushd SCCNV && unzip sccnv_example_v1.0.2.zip && bedtools makewindows -g resource/hg19.chrlength.txt -w 500000 > resource/hg19.bin500000.bed && popd

cp SCCNV/resource/hg19.mappability_gc.bin500000.bed SCCNV/resource/hg19.mappability_gc.bin500000.bed.old
for chr in chr $(seq 1 22) X Y; do cat SCCNV/resource/hg19.mappability_gc.bin500000.bed.old | grep -P ^$chr'\t' ; done | awk '{print "chr" $0 }' | sed 's/chrchr/chr/g' > SCCNV/resource/hg19.mappability_gc.bin500000.bed
cp SCCNV/resource/hg19.chrlength.txt SCCNV/resource/hg19.chrlength.txt.old
for chr in $(seq 1 22) X Y; do cat SCCNV/resource/hg19.chrlength.txt.old | grep -P ^$chr'\t' ; done | awk '{print "chr" $0 }' | sed 's/chrchr/chr/g' > SCCNV/resource/hg19.chrlength.txt
bedtools makewindows -g SCCNV/resource/hg19.chrlength.txt -w 500000 > SCCNV/resource/hg19.bin500000.bed

# https://github.com/deepomicslab/SeCNV
git clone https://github.com/deepomicslab/SeCNV.git
sed 's;(np.int);(int);g' SeCNV/Scripts/call_cn.py -i
pushd SeCNV/Scripts
wget -c https://zenodo.org/records/14407911/files/hg19_mappability.bigWig https://zenodo.org/records/14407911/files/hg38_mappability.bigWig
popd

#conda create --name chisel chisel
#mamba env create -f ../env/chisel.freeze.env_export.yml
conda create --yes --name chisel --file ../env/chisel.requirements.list_e_no_pypi.txt
#sed "s;profile = 'HS20' if abs(100 - read) < abs(125 - read) else 'HS25';profile = 'HS25' # 'HS20' if abs(100 - read) < abs(125 - read) else 'HS25';g" -i ${CONDA_PREFIX}/../chisel/lib/python2.7/site-packages/chisel/bin/chisel_nonormal.py
sed "s; -ss {} ;;g" -i ${CONDA_PREFIX}/../chisel/lib/python*/site-packages/chisel/bin/chisel_nonormal.py
sed "s;, profile, ;, ;g" -i ${CONDA_PREFIX}/../chisel/lib/python*/site-packages/chisel/bin/chisel_nonormal.py

# SCYN
for py in $(ls ${CONDA_PREFIX}/lib/python3.*/site-packages/scyn/scyn.py); do cp scyn/scyn.py $py ; done

