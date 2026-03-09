from itertools import zip_longest
import sys

umi_len=12
bc_len=16

umi_fastq = sys.argv[2]
bc_fastq = sys.argv[3]

# Line by line with zip_longest
def parse_zip_longest():
    fastq_iterator = (l.rstrip() for l in sys.stdin)
    for record in zip_longest(*[fastq_iterator] * 4):
        yield record

umi_file=open(umi_fastq, "w")
bc_file=open(bc_fastq, "w")

for read in parse_zip_longest():
	# First 16bp = visium bc, last 10bp = UMI
	# Write to bc file
	bc_file.write(f'{read[0]}\n')
	bc_file.write(f'{read[1][0:bc_len]}\n')
	bc_file.write(f'{read[2]}\n')
	bc_file.write(f'{read[3][0:bc_len]}\n')

	# Write to umi file
	umi_file.write(f'{read[0]}\n')
	umi_file.write(f'{read[1][bc_len:bc_len+umi_len]}\n')
	umi_file.write(f'{read[2]}\n')
	umi_file.write(f'{read[3][bc_len:bc_len+umi_len]}\n')