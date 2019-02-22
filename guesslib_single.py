import subprocess
import csv
import argparse
from guesslib_pair import run_blat
from guesslib_pair import find_good_blat_alignment_f1

###--------- GLOBAL VARIABLES ---------###

# Specify from what line of the user's transcripts
# that sequences should be aligned from.
START_FROM_SEQUENCE_NR = 3

# Specify from how many paired alignments data
# should be collected from.
NR_OF_PAIRS = 10

###--------- FUNCTIONS ---------###

def guesslib_single(ref, user_transcripts):
	'''
	Determines if the strand is forward, reverse or unstranded
	'''
	
	read_type = "N/A"
	read_start = START_FROM_SEQUENCE_NR * 4
	forward = 0
	reverse = 0

	for i in range(NR_OF_PAIRS):
		read_type = "N/A"

		(read_type,
			read_start) = find_and_determine_one_read(ref,
														user_transcripts,
														read_start)
		if (read_type == "F"):
			forward += 1
			print("Added to forward")
		elif (read_type == "R"):
			print("Added to reverse")
			reverse += 1

	# Determination of library type.
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

	print(f'The lib type is {lib_type}. Nr of forwards: {forward}. Nr of reverse: {reverse}.')

	return lib_type, forward, reverse

def find_and_determine_one_read(ref, user_transcripts, read_start):
	
	lib_type = "N/A"
	nr_of_results = 0

	while (nr_of_results != 1 or read_type == "N/A"):
		read_type = "N/A"

		(seq_id, seq_start, seq_end,
			nr_of_results, ref_transcript_id_R1,
			read_start) = find_good_blat_alignment_f1(ref,
                                        user_transcripts,
	                                            read_start)

		if (nr_of_results == 1):
			read_type = orientation_analysis(seq_start, seq_end)

	print(f'The read start and end is',
		f'({seq_start}, {seq_end})')

	return (read_type, read_start)

def orientation_analysis(seq_start, seq_end):
    '''
    This function determines the orientation of one read.
    '''

    read_type = "N/A"

    if (seq_start < seq_end):
    	read_type = "F"
    elif (seq_start > seq_end):
    	read_type = "R"

    return read_type

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
