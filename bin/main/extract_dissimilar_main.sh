#########################################################################################################################
# name: extract_dissimilar_main
# purpose: main script driving extract_dissimilar_transcripts
# input: $1 : reference transcriptome
#        $2 : a directory containing (and only containing) the transcriptomes to be added
#             (If you store these directories separately on your local machine, create a
#              directory with symbollic links to the appropriate transcriptomes)
# output: a concatenated transcriptome consisting of the reference transcriptome +
#         the most dissimilar sequences from each of the other transcriptomes
# written by: chase mateusiak, chase.mateusiak@gmail.com, chase.mateusiak@ucsf.edu
# date: 20190814
#
# Dependencies: blast+ > ___v, pandas etc from .py scripts(s), python3
#
# Note:
#########################################################################################################################

# main method takes input from cmd line, executes methods(camelCase, see methods section below)/scripts in repository bin
main($1, $2){
  # rename cmd line input for clarity
  local ref_transcriptome=$1
  local extra_transcriptome_dir=$2

  # get basenames. extra_transcriptome_basenames is a list.
  local ref_transcriptome_basename=$(getBasenames $extra_transcriptome_dir)
  local extra_transcriptome_basenames=$(getBasenames $extra_transcriptome_dir)

  # ask user to confirm/provide extra_transcriptome_names
  local ref_transcriptome_basename=$(userConfirmNames $ref_transcriptome_basename)
  local extra_transcriptome_basenames=$(userConfirmNames $extra_transcriptome_basenames)

  # create extract_dissimilar_transcripts results directory
  createDirectories $ref_transcriptome_basename $extra_transcriptome_basenames

  # make database from ref_transcriptome
  makeBlastDB $ref_transcriptome $ref_transcriptome_basename

  # blast ref_transcriptome x each extra_transcriptome
  againstRef_blastOutfmt6 $ref_transcriptome $extra_transcriptome_dir $extra_transcriptome_basenames

  # extract unique transcript names of each extra_transcriptome dir, deposit in appropriate subdir
  extractUniqTranscriptIDs /path/to/dir

  # create .fasta of most dissimilar transcripts in each transcriptome
  # transcripts not in the blastn ref_transcriptome x extra_transcriptome results,
  # which are transcripts which have no similarity)
  filterFasta /path/to/dir

  # create concat transcriptome .fasta MINUS GIVEN TRANSCRIPTOME of combination set {input transcriptomes}
  concatTranscriptomeMinus /path/to/dir

  # make blast database from each concat_minus, deposit in db appropriate subdir
  makeBlastDB

  # blast each filtered_extra_transcriptome against appropriate concat_transcriptome_minus_db
  againstFilteredTranscriptome_blastnOutfmt6

  # extract unique transcript IDs from filtered_transcriptomes, deposit in appropriate subdir
  extractUniqTranscriptIDs

  # filter the filtered_fasta
  filterFasta

  # make new fasta of the most dissimilar sequences to both reference as well as between extra filtered_transcriptomes
  concatTranscriptomes

} # end main

### step functions ###

getBasenames(){

} # end getBasenames()

userConfirmNames(){
  # If user renames, create list of tuples (basename, username)

} # end userConfirmNames()

createDirectories(){
  # incl subdir db for each extra_transcriptome

} # end createDirectories()

makeBlastDB(){
  # put make_blast_db.sh $ref_transcriptome $ref_transcriptome_basename here
} # end makeBlastDB

againstRef_blastOutfmt6(){

} # end againstRef_blastOutfmt6()

extractUniqTranscriptIDs(){

} # end extractUniqTranscriptIDs()

filterFasta(){
  create_dissimilar_transcriptome.py > ${basename}_blah.fasta
} # end createFilterFasta()

concatTranscriptomeMinus(){

} # end createConcatMinus()

### helper functions ###
againstFilteredTranscriptome_blastnOutfmt6(){
  # put blastn_output6.sh $ref_transcriptome $extra_transcriptome_dir $extra_transcriptome_basenames here
} # end makeBlastDB()

concatTranscriptomes(){

} # end concatTranscriptomes
