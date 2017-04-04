from Bio import SeqIO
input_seq_iterator = SeqIO.parse("all_segment.fa", "fasta")
virus_seq_iterator = (record for record in input_seq_iterator \
                      if "=1:/10239:0/" in record.description)


SeqIO.write(virus_seq_iterator, "virus_segment.fa", "fasta")

input_seq_iterator = SeqIO.parse("all_segment.fa", "fasta")

nonvirus_seq_iterator = (record for record in input_seq_iterator \
                      if "=1:/10239:0/" not in record.description)


SeqIO.write(nonvirus_seq_iterator, "nonvirus_segment.fa", "fasta")
