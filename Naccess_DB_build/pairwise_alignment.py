from Bio import pairwise2
from Bio import SeqIO
import glob
import pandas as pd

with open("./template_sequences/template_names.log", 'w') as f:
	for i, filename in enumerate(glob.glob('./template_sequences/allele/*.pdb')):
		seq2_alpha = SeqIO.read(filename, "fasta")
		f.write(seq2_alpha.id + "\n")
		print(seq2_alpha.id)
		best_record_list = []
		best_score = 0
		with open("MHC_seq.fasta") as handle:
			for record in SeqIO.parse(handle, "fasta"):
				alignments = pairwise2.align.globalxx(seq2_alpha.seq, record.seq)
				if(alignments[0].score == best_score):
					best_record_list.append(record.id)
				if(alignments[0].score > best_score):
					best_score = alignments[0].score
					best_record_list = []
					best_record_list.append(record.id)
			f.write(best_record_list[0] + "\n")
			print(best_record_list[0])
f.close()