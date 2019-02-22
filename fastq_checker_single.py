import argparse
import os
import subprocess

###--- GLOBAL VARIABLES ---###

# Specify the least number of sequences in the out-file
NR_OF_SEQUENCES_THRESHOLD = 1000

# Specify the quality threshold you want these sequences to have
# If you don't get NR_OF_SEQUENCES seqs of this quality, the quality will
# lower in order to accomodate that
QUALITY_THRESHOLD_START = 30

###--- FUNCTIONS ---###

def proper_fastq_format(fastq_file):
	'''
	This function returns True if the fastq file is properly formatted.
	'''

	fastq_format = False
	fastq_decision = 0
	counter = 0

	with open(fastq_file) as f_in:
		for line in f_in:
			if (counter % 4 == 0): # Checking first line
				if (line.startswith('@')):
					fastq_decision += 1
			if (counter % 4 == 1): # Checking second line
				second_line = line
				if all(c in 'ATGCNatgcn\n' for c in list(line)):
					fastq_decision += 1
			if (counter % 4 == 2): # Checking third line
				if (line.startswith('+')):
					fastq_decision += 1
			if (counter % 4 == 3): # Checking fourth line
				if (len(line) == len(second_line)):
					fastq_decision += 1

			counter += 1

	if (fastq_decision == counter):
		fastq_format = True

	return fastq_format


def calculate_phred_quality(phred_string):
	'''
	This function calculates the average Q-phred based on an ascii-string,
	using the phred ascii-33 system.
	'''

	phred_char_list = list(phred_string)
	phred_quality = 0

	for i in range(len(phred_char_list)):
		phred_quality = phred_quality + ord(phred_char_list[i])-33

	phred_quality = phred_quality/len(phred_char_list)
	
	return phred_quality


def remove_low_quality_reads(fastq_file, quality_threshold):
	'''
	This function removes low quality reads as specified by
	an average Q-threshold.
	'''

	nr_of_sequences = 0
	q_fastq_file = "Q" + str(int(quality_threshold)) + "_" + fastq_file

	f = open(q_fastq_file, "w+")
	f.close()

	with open(fastq_file) as f_in:
		counter = 0
		for line in f_in:
			if (counter % 4 == 0):
				first_line = line
			if (counter % 4 == 1):
				second_line = line
			if (counter % 4 == 2):
				third_line = line
			if (counter % 4 == 3):
				fourth_line = line
				phred_quality = calculate_phred_quality(fourth_line)

				if (phred_quality > quality_threshold):
					nr_of_sequences += 1
					with open(q_fastq_file, "a") as f_out: 
						f_out.write(first_line)
						f_out.write(second_line)
						f_out.write(third_line)
						f_out.write(fourth_line)

			counter += 1

	return q_fastq_file, nr_of_sequences, quality_threshold


def remove_low_quality_reads_without_q(fastq_file):
	'''
	This function removes low quality reads without a specified Q-theshold.
	Instead, it will remove the low quality reads below QUALITY_THRESHOLD_START.
	If the number of sequences are below NR_OF_SEQUENCES_THRESHOLD, it will
	iterate with a lower threshold.
	'''

	# Initializing values
	quality_threshold = QUALITY_THRESHOLD_START
	nr_of_sequences_threshold = NR_OF_SEQUENCES_THRESHOLD
	nr_of_sequences = 0
	q_fastq_file = "N/A"
	proper_format = False

	if proper_fastq_format(fastq_file):
		proper_format = True

		(q_fastq_file,
			nr_of_sequences,
			quality_threshold) = remove_low_quality_reads(fastq_file,
														quality_threshold)

		while (nr_of_sequences < nr_of_sequences_threshold):
			if os.path.exists(q_fastq_file):
				subprocess.run(['rm', q_fastq_file]) # Remove the file generated
														# in previous iteration
			quality_threshold = quality_threshold - 1 # Lower Q-treshold

			(q_fastq_file,	# Get new files file-pair
				nr_of_sequences,
				quality_threshold) = remove_low_quality_reads(fastq_file,
															quality_threshold)

	return (q_fastq_file, nr_of_sequences,
			quality_threshold, proper_format)


###--- MAIN ---###

def main():

	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--file", help="FASTQ file")
	args = parser.parse_args()

	if not len(vars(args)) == 1:
		print("Specify the file and quality threshold")

	else:
			(q_fastq_file, nr_of_sequences,
				quality_threshold,
				proper_format) = remove_low_quality_reads_without_q(args.file)
			
			print(f'The input file was in proper fastq format: {proper_format}\n',
				f'The subsetted file is called: {q_fastq_file}',
				f'It contains {nr_of_sequences} nr of sequences',
				f'with quality above Q={quality_threshold}.')


if __name__ == "__main__":
	main()
