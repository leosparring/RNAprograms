from fastq_checker_single import remove_low_quality_reads
from fastq_checker_single import proper_fastq_format
import subprocess
import argparse
import os

###--- GLOBAL VARIABLES ---###

# Specify the least number of sequences in the out-files
NR_OF_SEQUENCES_THRESHOLD = 1000

# Specify the quality threshold you want these sequences to have
# If you don't get NR_OF_SEQUENCES seqs of this quality, the quality threshold will
# lower in order to accomodate that
QUALITY_THRESHOLD_START = 40

###--- FUNCTIONS ---###

def subset_to_pairs_without_q(fastq_file_1, fastq_file_2):
	'''
	This function subsets the two fastq files and returns at least
	the nr of sequences specified as the global variable. It starts with
	the quality threshold specified in the second global variable.
	The function returns the nr of sequences in the subsetted files, the
	quality threshold used, and whether the subsetting was succesful
	or not.
	'''
	
	quality_threshold = QUALITY_THRESHOLD_START
	nr_of_sequences_threshold = NR_OF_SEQUENCES_THRESHOLD

	(subset_successful, 
		nr_of_sequences,
		new_fastq_file_1,
		new_fastq_file_2) = subset_to_pairs(fastq_file_1,
											fastq_file_2,
											quality_threshold)

	while (subset_successful and nr_of_sequences < nr_of_sequences_threshold):
		if os.path.exists(new_fastq_file_1):
			subprocess.run(['rm', new_fastq_file_1])
			subprocess.run(['rm', new_fastq_file_2])

		quality_threshold = quality_threshold - 1
		
		(subset_successful, 
		nr_of_sequences,
		new_fastq_file_1,
		new_fastq_file_2) = subset_to_pairs(fastq_file_1,
											fastq_file_2,
											quality_threshold)

	return (subset_successful, nr_of_sequences,
			quality_threshold, new_fastq_file_1, new_fastq_file_2)


def subset_to_pairs(fastq_file_1, fastq_file_2, quality_threshold):
	'''
	This function subsets the two fastq files to the specified threshold,
	then it removes the transcripts that are only in one of the two files.
	In this way there will still be a complete set of pairs after subsetting.
	'''

	proper_format = False
	nr_of_sequences = 0
	pairchecked_file_1 = "N/A"
	pairchecked_file_2 = "N/A"

	if (proper_fastq_format(fastq_file_1)
		and proper_fastq_format(fastq_file_2)):
		proper_format = True

		q_fastq_file_1 = remove_low_quality_reads(fastq_file_1,
													quality_threshold)[0]
		q_fastq_file_2 = remove_low_quality_reads(fastq_file_2,
													quality_threshold)[0]
		
		# Naming outfiles
		pairchecked_file_1 = 'paircheck_' + q_fastq_file_1
		pairchecked_file_2 = 'paircheck_' + q_fastq_file_2

		# Initiating lists of sequences
		seq_ids_1 = []
		seq_ids_2 = []

		# Get a list of sequence ID's in q_fastq_file_1
		with open(q_fastq_file_1) as fin_1:
			counter = 0
			for line in fin_1:
				if (counter % 4 == 0):
					seq_id = line.split(" ")[0]
					seq_ids_1.append(seq_id)
				counter += 1

		# Get a list of sequence ID's in q_fastq_file_2.
		# Check if the sequences in q_fastq_file_2 are in the list seq_ids_1
		# if they are, write them to pairchecked_file_2
		# Also, determine the number of sequences in the final files.
		with open(pairchecked_file_2, 'w') as f_out:
			with open(q_fastq_file_2) as fin_2:
				counter = 0
				for line in fin_2:
					if (counter % 4 == 0):
						first_line = line
						seq_id = first_line.split(" ")[0]
						seq_ids_2.append(seq_id)
						if seq_id in seq_ids_1:
							add_to_new = True
							nr_of_sequences += 1
							f_out.write(first_line)
						else:
							add_to_new = False
					if (add_to_new):
						if (counter % 4 == 1):
							f_out.write(line)
						if (counter % 4 == 2):
							f_out.write(line)
						if (counter % 4 == 3):
							f_out.write(line)
							add_to_new = False	
					counter += 1

		# Check if the sequences in q_fastq_file_2 are in the list seq_ids_2.
		# If they are, write them to pairchecked_file_1.
		with open(pairchecked_file_1, 'w') as f_out:
			with open(q_fastq_file_1) as fin_1:
				counter = 0
				for line in fin_1:
					if (counter % 4 == 0):
						first_line = line
						seq_id = first_line.split(" ")[0]
						if seq_id in seq_ids_2:
							add_to_new = True
							f_out.write(first_line)
						else:
							add_to_new = False
					if (add_to_new):
						if (counter % 4 == 1):
							f_out.write(line)
						if (counter % 4 == 2):
							f_out.write(line)
						if (counter % 4 == 3):
							f_out.write(line)
							add_to_new = False	
					counter += 1

		# Remove the files used as input to this function.
		subprocess.run(['rm', q_fastq_file_1])
		subprocess.run(['rm', q_fastq_file_2])

	return (proper_format, nr_of_sequences,
			pairchecked_file_1, pairchecked_file_2)


###--- MAIN ---###

def main():

	parser = argparse.ArgumentParser()
	parser.add_argument("-f1", "--file_1", required=True, help="FASTQ file one")
	parser.add_argument("-f2", "--file_2", required=True, help="corresponding FASTQ of paired sequences")
	args = parser.parse_args()

	
	(proper_format,
		nr_of_sequences,
		quality,
		out_file_1,
		out_file_2) = subset_to_pairs_without_q(args.file_1,
												args.file_2)

	print(f'The two input files are properly formatted: {proper_format}\n',
		f'The number of sequences in your subsetted files are: {nr_of_sequences}\n',
		f'All of these sequences have a quality higher than Q={str(quality)}\n',
		f'The file names are {out_file_1} and {out_file_2}')

if __name__ == "__main__":
	main()