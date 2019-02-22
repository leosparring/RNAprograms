import subprocess
import csv
import argparse

###--------- GLOBAL VARIABLES ---------###

# Specify from what line of the user's transcripts
# that sequences should be aligned from.
START_FROM_SEQUENCE_NR = 0

# Specify from how many paired alignments data
# should be collected from.
NR_OF_PAIRS = 50

###--------- FUNCTIONS ---------###

def guesslib_pair(ref, user_transcripts_f1, user_transcripts_f2):
    '''
    Determines if the strand is forward, reverse or unstranded, based
    on the nr of pairs that indicate pair-types.
    '''

    lib_type = 'N/A'
    read_start = START_FROM_SEQUENCE_NR * 4
    o_and_f = 0
    o_and_r = 0
    i_and_f = 0
    i_and_r = 0
    invalid_orientation_strand = 0

    # Collecting pairs
    for i in range(NR_OF_PAIRS):
        pair_type = "N/A"

        (pair_type, read_start) = find_a_significant_pair(ref,
                                            user_transcripts_f1,
                                            user_transcripts_f2,
                                            read_start)

        if (pair_type == "OF"):
            o_and_f += 1
            print("Added to OF.\n")
        elif (pair_type == "OR" ):
            o_and_r += 1
            print("Added to OR.\n")
        elif (pair_type == "IF"):
            i_and_f += 1
            print("Added to IF.\n")
        elif (pair_type == "IR"):
            i_and_r += 1
            print("Added to IR.\n")
        else:
            invalid_orientation_strand += 1
            print("Added to none of them.\n")

    # Determination of library type.
    if ((i_and_r + i_and_f) > (o_and_r + o_and_f)):
        if (i_and_r == 0):
            ratio_IF_over_IR = 9
        else:
            ratio_IF_over_IR = (i_and_f/i_and_r)
        if (ratio_IF_over_IR > 8):
            lib_type = 'ISF'
        elif (ratio_IF_over_IR < (1/8)):
            lib_type = 'ISR'
        else:
            lib_type = 'IU'
    if ((i_and_r+i_and_f) < (o_and_r+o_and_f)):
        if (o_and_r == 0):
            ratio_OF_over_OR = 9
        else:
            ratio_OF_over_OR = (o_and_f/o_and_r)
        if (ratio_OF_over_OR > 8):
            lib_type = 'OSF'
        elif (ratio_OF_over_OR < (1/8)):
            lib_type = 'OSR'
        else:
            lib_type = 'UO'

    print(f'\nOF: {o_and_f}, OR: {o_and_r}',
            f'IF: {i_and_f}, IR: {i_and_r}, neither: {invalid_orientation_strand}.',
            f'The most likely lib type is {lib_type}.')

    return (lib_type, o_and_f, o_and_r, i_and_f, i_and_r)


def pair_analysis(seq_start_R1, seq_end_R1, seq_start_R2, seq_end_R2):
    '''
    This function determines inward or outward and reverse
    or forward orientation based on the input of two reads.
    '''

    pair_type = "N/A"

    distance_r2_r1 = (((seq_start_R2 + seq_end_R2)/2) -
                    ((seq_start_R1 + seq_end_R1)/2))

    if (seq_start_R1 < seq_end_R1 and seq_start_R2 > seq_end_R2):
        if (distance_r2_r1 > 0 and distance_r2_r1 < 1000):
            pair_type = "IF" #
        elif (distance_r2_r1 < -750):
            pair_type = "OF"
    elif (seq_start_R1 > seq_end_R1 and seq_start_R2 < seq_end_R2):
        if (distance_r2_r1 < 0 and distance_r2_r1 > -1000):
            pair_type = "IR"
        elif (distance_r2_r1 > 750):
            pair_type = "OR"

    print(f'(Avg of R2) - (Avg of R1) is: {distance_r2_r1}\n')

    return pair_type


def find_a_significant_pair(ref, user_transcripts_f1, user_transcripts_f2, read_start):
    '''
    This function finds a significant alignment in transcript file 1,
    and then finds the corresponding sequence's alignment in file 2.
    If they are not significant, the first fastq file is searched further.
    '''

    nr_of_results_R1 = 0
    nr_of_results_R2 = 0
    pair_type = "N/A"

    while ((nr_of_results_R1 != 1 or nr_of_results_R2 != 1)
            or (pair_type == "N/A")):

        nr_of_results_R1 = 0
        nr_of_results_R2 = 0 # Resetting the hits

        (seq_id, seq_start_R1, seq_end_R1,
            nr_of_results_R1, ref_transcript_id_R1,
            read_start) = find_good_blat_alignment_f1(ref,
                                                user_transcripts_f1,
                                                read_start)

        (seq_start_R2, seq_end_R2, nr_of_results_R2,
            ref_transcript_id_R2) = blat_corresponding_seq_in_f2(ref,
                                                user_transcripts_f2,
                                                seq_id)

        if (nr_of_results_R1 == 1 and nr_of_results_R2 == 1
            and ref_transcript_id_R1 == ref_transcript_id_R2):
            pair_type = pair_analysis(seq_start_R1,
                                        seq_end_R1,
                                        seq_start_R2,
                                        seq_end_R2)

    print(f'The R1 start and end is',
            f'({seq_start_R1}, {seq_end_R1}). ',
            f'The R2 start and end is ({seq_start_R2},{seq_end_R2}).')

    return (pair_type, read_start)


def blat_corresponding_seq_in_f2(ref, user_transcripts, seq_id_R1):
    '''
    This function finds the corresponding seq to the id of the input in the
    reference fastq file specified as ref. It then performs a blat and returns
    the start, end and e-value of this hit.
    '''

    seq_id_R2 = "N/A"
    iterations = 0
    ref_transcript_id = "N/A"
    tmp_in_2 = "tmp/tmp_blat_inpt_2_" + user_transcripts # Naming temp fasta file

    with open(user_transcripts) as f_in:
        first_line = '>' + f_in.readline()[1:]
        seq_id_R2 = first_line.split(" ")[0]

        while (seq_id_R2 != seq_id_R1 and iterations < 1000):
            next(f_in)
            next(f_in)
            next(f_in)
            first_line = '>' + f_in.readline()[1:]
            seq_id_R2 = first_line.split(" ")[0]
            iterations += 1

        with open(tmp_in_2, "w+") as f_out:
            f_out.write(first_line)
            f_out.write(f_in.readline())

    print("Now blatting the corresponding sequence in file 2.")

    (seq_start, seq_end, nr_of_results,
        ref_transcript_id) = run_blat(ref, tmp_in_2)

    subprocess.run(['rm', tmp_in_2]) # Removing the tmp file

    return (seq_start, seq_end,
            nr_of_results, ref_transcript_id)

def find_good_blat_alignment_f1(ref, user_transcripts, read_start):
    '''
    This function reads the user's first transcript file line by line,
    starting from line read_start, and extracts sequences. Then it runs
    a blat search with the sequence. If it meets the conditions, it returns
    the seq id, start and end, and the significance of the result.
    '''

    nr_of_results = 0
    iterations = 0

    tmp_in_1 = "tmp/tmp_blat_inpt_1_" + user_transcripts #naming temp fasta file

    with open(user_transcripts) as f_in:
        for i in range(read_start):
            next(f_in)
        while (nr_of_results != 1):
            with open(tmp_in_1, "w+") as f_out: # Creating FASTA tmp
                first_line = '>' + f_in.readline()[1:]
                seq_id = first_line.split(" ")[0]
                print(seq_id)
                f_out.write(first_line)  # Writes the header into f_out.
                f_out.write(f_in.readline()) # Writes the sequence into f_out.
            next(f_in)
            next(f_in)
            (seq_start, seq_end, nr_of_results,
                ref_transcript_id) = run_blat(ref, tmp_in_1)

            iterations += 1

    read_start = read_start + iterations * 4

    subprocess.run(['rm', tmp_in_1]) # Removing the tmp file

    return (seq_id, seq_start, seq_end,
            nr_of_results, ref_transcript_id, read_start)


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

    subprocess.run(['rm', tmp_rslt]) # Removing the tmp file

    return (int(seq_start), int(seq_end),
            nr_of_results, ref_transcript_id)

###--------- MAIN ---------###

def main():
	
	parser = argparse.ArgumentParser()
	parser.add_argument("-r", "--reference", required=True, help="reference sequence")
	parser.add_argument("-f1", "--file_1", required=True, help="FASTQ file one")
	parser.add_argument("-f2", "--file_2", required=True, help="corresponding FASTQ of paired sequences")
	args = parser.parse_args()

	if not len(vars(args)) == 3:
		print("Specify the file 2 fastq files and reference library")

	else:
		guesslib_pair(args.reference, args.file_1, args.file_2)

if __name__ == "__main__":
	main()
