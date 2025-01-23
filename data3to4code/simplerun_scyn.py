import sys
import scyn

bam_dir = sys.argv[1]
output_dir = sys.argv[2]
csv = sys.argv[3]

# create SCYN object
scyn_operator = scyn.SCYN(reg='*.bam$', seq='paired-end', bin_len = 1000)

# call cnv
# bam_dir is the input bam directory and output_dir is the output directory
scyn_operator.call(bam_dir, output_dir)

# store cnv matrix to a csv file
scyn_operator.cnv.to_csv(csv)
