# RNA-seq Programs

These programs are used to subset fastq data and subsequently determine their library type by blatting
the sequences to a transcriptome. They are intended to be used as two
separate pipelines. They were developed by students in Uppsala university
bioinformatics program. The purpose of the development was for use for the online library type
determinator 'guesslib'.

Note: These programs have only been tested with linux. In order to
use these programs you need a blat executable installed in the same folder.
Download blat from here:
http://hgdownload.soe.ucsc.edu/admin/exe/
You will also need to create a subdirectory "tmp" and "reference_sequences". Place your
transcriptome sequences in the latter of these.

# Usage
## 1
subset your raw fastq file(s) to the first 10000 lines,
e.g. in linux use the command head -n 10000 <in_filename> > <out_filename>.

## 2
Perform a quality check on your fastq file(s) by using the fastq_checker_single.py for
single end files, and fastq_checker_pair.py for paired end, split fastq files.
In the terminal:
python3 fastq_checker_single.py -f <subsetted_fastq>
python3 fastq_checker_single.py -f1 <subsetted_f1.fastq> -f2 <subsetted_f2.fastq>

## 3
determine the library type by running the guesslib_pair.py or guesslib_single.py programs.
In the terminal:
python3 guesslib_single.py -f <subsetted_fastq> -r <reference_transcriptome>
python3 guesslib_pair.py -f1 <subsetted_f1.fastq> -f2 <subsetted_f2.fastq> -r <reference_transcriptome>

