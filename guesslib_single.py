import subprocess
import csv
import argparse

###--------- GLOBAL VARIABLES ---------###

# Specify from what line of the user's transcripts
# that sequences should be aligned from.
START_FROM_SEQUENCE_NR = 0

# Specify from how many paired alignments data
# should be collected from.
NR_OF_PAIRS = 2

# Specify the maximum nr of alignments to perform before giving up on determining
# the library type
MAX_BLATS = 200

###--------- FUNCTIONS ---------###

def guesslib_single(ref, user_transcripts):
	'''
	Finds NR_OF_PAIRS pairs and
	determines if the library is forward, reverse or unstranded
	'''

	# Initializing variables
	collected_reads = 0
	read_start = START_FROM_SEQUENCE_NR * 4
	forward = 0
	reverse = 0
	# Naming tmp fasta file according to input file name
	tmp_fasta = "tmp/tmp_blat_inpt_" + user_transcripts
	seqs_searched = 0
	lib_type = "N/A"
	succesful_lib_determination = False

	with open (user_transcripts) as f_in:
		for i in range(read_start): # Skipping to START_FROM_SEQUENCE_NR
			next(f_in)
		while (collected_reads < NR_OF_PAIRS and seqs_searched < MAX_BLATS):
			read_type = "N/A"
			nr_of_results = 0

			write_tmp_fasta(tmp_fasta, f_in)
			(seq_start, seq_end, nr_of_results,
				ref_transcript_id) = run_blat(ref, tmp_fasta)
			seqs_searched += 1
	
			if (nr_of_results == 1):
				read_type = orientation_analysis(seq_start, seq_end)
				if (read_type == "F"):
					forward += 1
				elif (read_type == "R"):
					reverse += 1
				collected_reads += 1
				print("Collected reads: " + str(collected_reads) + "\n")

	subprocess.run(['rm', tmp_fasta]) # Removing the tmp FASTA file

	if (seqs_searched == MAX_BLATS):
		print(f'Could not determine library type\n'
			f'Tried to align {MAX_BLATS} sequences'
			f' of which {collected_reads} aligned according to the criteria')
	else:
		lib_type = determine_lib_type(forward, reverse)
		succesful_lib_determination = True
	
	return lib_type, forward, reverse, succesful_lib_determination

def determine_lib_type(forward, reverse):
	'''
	This function determines the library type, based on
	the number of forward reads in relation to the
	number of reverse reads.
	'''

	if (reverse == 0):
		ratio_F_over_R = 9
	else:
		ratio_F_over_R = (forward/reverse)
	if (ratio_F_over_R > 8):
	    lib_type = 'SF'
	elif (ratio_F_over_R < (1/8)):
	    lib_type = 'SR'
	else:
	    lib_type = 'U'

	print(f'The lib type is {lib_type}.'
			f' Nr of forwards: {forward}. Nr of reverse: {reverse}.')

	return (lib_type)

def orientation_analysis(seq_start, seq_end):
	'''
	This function determines the orientation of one read.
	'''

	read_type = "N/A"

	print(f'The read start and end is ({seq_start}, {seq_end})')
	if (seq_start < seq_end):
		read_type = "F"
		print("The read type is F")
	elif (seq_start > seq_end):
		read_type = "R"
		print("The read type is R")

	return read_type

def write_tmp_fasta(tmp_fasta, fastq):
	'''
	Writes a temporary fasta file, taking in a fastq-file
	'''

	with open(tmp_fasta, "w+") as f_out: # Creating FASTA tmp
		first_line = '>' + fastq.readline()[1:]
		seq_id = first_line.split(" ")[0]
		print(seq_id)
		f_out.write(first_line)  # Writes the header into f_out.
		f_out.write(fastq.readline()) # Writes the sequence into f_out.
		next(fastq) # Skipping past the +, and phred lines.
		next(fastq)

def run_blat(ref, seq):
    ''' This function calls blat with a single sequence seq.fa as
        input against the reference ref.fa and returns seq_identity,
        seq_start, seq_end'''

    nr_of_results = 0
    seq_identity = "0"
    seq_start = "0"
    seq_end = "0"
    second_hit_e_value = "0"
    e_value = "0"
    ref_transcript_id = "N/A"
    tmp_rslt = seq + "_rslt"
    ref = 'reference_sequences/' + ref

    subprocess.run(['./blat', ref, seq, '-out=blast8', tmp_rslt])
    for line in open(tmp_rslt):
        nr_of_results += 1
    if (nr_of_results > 0):
        with open(tmp_rslt, 'r') as f_in:
            first_hit = next(csv.reader(f_in, delimiter='\t'))
            ref_transcript_id = first_hit[1]
            seq_identity = first_hit[2]
            seq_start = first_hit[8]
            seq_end = first_hit[9]
            e_value = first_hit[10]
            if (nr_of_results > 1):
                second_hit = next(csv.reader(f_in, delimiter='\t'))
                second_hit_e_value = second_hit[10]

    print(f'the identity is: {seq_identity}')
    print(f'the E value of the first hit is: {e_value}')
    print(f'the seq start is: {seq_start}')
    print(f'the seq end is: {seq_end}')
    print(f"the second hit's E value is: {second_hit_e_value}")
    print(f'the reference transcript id is {ref_transcript_id}')
    print(f'the number of results for this blat is: {nr_of_results}\n')

    subprocess.run(['rm', tmp_rslt]) # Removing the tmp blast-rslt file

    return (int(seq_start), int(seq_end),
            nr_of_results, ref_transcript_id)

###--- MAIN ---###

def main():
	
	parser = argparse.ArgumentParser()
	parser.add_argument("-r", "--reference", required=True, help="reference sequence")
	parser.add_argument("-f", "--user_transcripts", required=True, help="FASTQ file one")
	args = parser.parse_args()

	if not len(vars(args)) == 2:
		print("Specify the fastq file and reference library")

	else:
		guesslib_single(args.reference, args.user_transcripts)

if __name__ == "__main__":
	main()
