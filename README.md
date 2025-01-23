
step$(N)data is taken as input to step$(N)soft to output step$(N+1)data using step$(N)temp as the temporary directory

for a set of cells derived from the same donor (e.g., human subject):
    data1: sra-raw-fastq
      run data1to2bin
    data2: ref-aligned-bam with two germline haplotypes
      run data2to3bin (oocyte-variation-1: merging (PB1, PB2, FPN) data; oocyte-variation-2: separating (PB1, PB2, FPN) data)
    data3: ref-aligned-bam with simulated copy numbers
      run data3to4bin
    data4: cnv-call-tsv
      run data4to5bin
    data5: cnv-eval-tsv

#data  after step5: figures-and-tables

data1
data\_2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_4_<tool>

1to2.2_<donor>.step<...>.sh
1to2.2_<donor>.step<...>.tmpdir/
1to2_2_<donor>/

2to3.2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>.sh
2to3.2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>.tmpdir
2to3_2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>/

2to4.2_<donor>_3_<sampleType>_<avgSpotLen>_4_<tool>.sh
2to4.2_<donor>_3_<sampleType>_<avgSpotLen>_4_<tool>.tmpdir
2to4_2_<donor>_3_<sampleType>_<avgSpotLen>_4_<tool>/

3to4.2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_4_<tool>.sh
3to4.2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_4_<tool>.tmpdir
3to4_2_<donor>_3_<sampleType>_<avgSpotLen>_<cellLine>_4_<tool>/

<donor>.data1to2.step${N}....sh
<donor>.data1to2.tmpdir/...
<donor>.data2/<acc>...
<donor>\_<sampleType>\_<avgSpotLen>.data2to4\_<tool>.sh
<donor>\_<sampleType>\_<avgSpotLen>.data2to4\_<tool>.tmpdir/...
<donor>\_<sampleType>\_<avgSpotLen>.data4\_<tool>/<acc>\_<tool>...
<donor>\_<sampleType>\_<avgSpotLen>\_<cellLine>.data2to3.sh
<donor>\_<sampleType>\_<avgSpotLen>\_<cellLine>.data2to3.tmpdir/...
<donor>\_<sampleType>\_<avgSpotLen>\_<cellLine>.data3/<acc_1>\_<acc_2>\_<cellLine>
<donor>\_<sampleType>\_<avgSpotLen>\_<cellLine>.data3to4\_<tool>.sh
<donor>\_<sampleType>\_<avgSpotLen>\_<cellLine>.data3to4\_<tool>.tmpdir
<donor>\_<sampleType>\_<avgSpotLen>\_<cellLine>.data4\_<tool>/<acc_1>\_<acc_2>\_<cellLine>\_<tool>

