#########################################################################################################################
# name: extract_dissimilar_main
# purpose: main script driving extract_dissimilar_transcripts
# input: Run script from location at which you wish to deposit the results directory
#        $1 : reference transcriptome
#        $2 : a directory containing (and only containing) the transcriptomes to be added (excluding reference transcriptome).
#             If you store these directories separately on your local machine, create a
#              directory with symbollic links to the appropriate transcriptomes
#        $3 : project name
# output: a concatenated transcriptome consisting of the reference transcriptome +
#         the most dissimilar sequences from each of the other transcriptomes
# written by: chase mateusiak, chase.mateusiak@gmail.com, chase.mateusiak@ucsf.edu
# date: 20190814, edit 20190815
#
# Dependencies: blast+ > ___v, pandas etc from .py scripts(s), python3
#
# Note: .fa must be standard -- pacbio names transcriptome needed to be manipulated to
#########################################################################################################################

main(){
  # rename cmd line input for clarity
  local reference_transcriptome=$1
  local input_transcriptomes=$2
  local project_name=$3

  # initialize array to hold transcript .fa to add to reference_transcriptome at end of process
  local dissimilar_set=()
  local transcripts_to_concat=()

  # create project directory at current location, cd, store path
  makeDirChangeDir extract_dissimilar_${proj_name}
  local project_dir=$(pwd)

  # make database from reference_transcriptome
  makeDirChangeDir reference_database
  makeBlastDB $reference_transcriptome
  # store path to reference db (/path/to/dir/basename_of_database_files.*)
  local reference_db_name=$(ls .)
  local reference_transcriptome_db=$(realpath $reference_db_name)/${reference_db_name}

  # return to project directory
  cd $project_dir

  # create dir new_transcriptomes which will hold data for each input_transcriptome, cd, store path
  mkdirCd input_transcriptomes
  local input_transcriptomes_dir=$(pwd)

  # loop through input_transcriptomes, blast transcript against ref, get unique names of transcripts found to be similar to reference,
  # use unique names to extract transcripts without any similarity to reference and create .fa of these
  for transcriptome in $input_transciptomes;
    do
      # create directory and subdirectories for each input_transcriptome with subdir blast_against_ref, cd into it
      mkdirCd $(basename $transcriptome .fa)
      mkdirCd blast_against_reference
      # extract transcripts without similarities to ref_transcriptome, output as .fa
      makeDissimilarFasta $reference_transcriptome_db $transcriptome
      # add .fa to dissimilar_set
      dissimilar_set+=$(realpath $(find . -name '*.fa'))
      # move out of individual input_transcriptome to input_transcriptome dir
      cd $input_transcriptomes_dir
    done

  # loop through each new_transcriptome dir, create concat_transcriptome minus the
  # individual new transcriptome, create new blast database from this concatMinus db
  # blastn the dissimilar .fa against concatMinus_db and then extract those dissimilar transcriptomes,
  # which will be the unique transcripts from the given new_transcriptome to be added to the ref_transcriptome
  # in the final concatenation

  for input_transcriptome_dir in $(ls $input_transcriptomes_dir);
    do
      cd $input_transcriptome_dir
      # make new directory to hold data related to within_dissimilar set comparison
      mkdirCd blast_against_dissimilar
      # current input_transcriptome dissimilar_to_reference .fa
      local dissimilar_against_reference_transcriptome=$(realpath $(find ../blast_against_reference -name '*.fa'))
      # input dissimilar .fa in input_transcriptome_dir and array of all dissimilar .fa, create concat .fa of all dissimilar_input_transcriptomes EXCEPT current input_transcriptome
      mkdirCd concat_minus_fa
      createDissimilarConcatMinusFasta $dissimilar_against_reference_transcriptome $dissimilar_set
      cd ..
      # make blast database from concatMinus .fa
      makeBlastDB $(find ./concat_minus -name '*.fa')
      # store path to concatMinus database
      local databaseMinus_name=$(find . -name '*_db')
      local databaseMinus_path=$(realpath $database_name/${database_name})
      # create .fa of transcripts which dissimilar to the transcripts which are also dissimilar to the reference transcriptome
      makeDissimilarFasta $databaseMinus_path $dissimilar_against_reference_transcriptome
      # add this .fa to array which stores final .fa to concat to reference_transcriptome
      transcripts_to_concat+=$(realpath new dissimilar .fa)
      cd $input_transcriptomes_dir
    done

  cd $project_dir
  createFinalConcat $ref_transcriptome $transcripts_t_concat
} # end main()

mkdirCd(){
  mkdir $1
  cd $1
} # end mkdirCd

makeBlastDB (){
  local fasta=$1
  local db_name=$(basename $fa_to_db .fa)
  mkdirCd $db_name
  makeblastdb -parse_seqids -dbtype nucl -in $fasta -out $db_name
  cd ..
} # end makeBlastDB

makeDissimilarFa(){

  local db=$1
  local query_fasta=$2
  local query_fasta_bn=$(basename $2 .fa)
  local output_tsv=${query_fasta_bn}.tsv

  # blast query_fasta against db
  blastn -outfmt 6 -num_threads 8 -db $db -query $query_fasta -out $output_tsv

  #extract unique names from output of above
  awk -v FS='\t' '{print $1}' $output_tsv | uniq > $(basename output_tsv .tsv)_unique.tsv

  # create .fa of transcripts from $query_fasta not found to have any similarity to transcripts in $db
  if [[ $query_fasta == p* ]]; then
    create_fa_blast_unmatched.py $query_fasta $(basename output_tsv .tsv)_unique.tsv
  else
    create_fa_blast_unmatched.py $query_fasta $(basename output_tsv .tsv)_unique.tsv
  fi
} # end makeDissimilarFa()

createDissimilarConcatMinusFasta(){
  local dissim_against_ref_fasta=$1
  local dissim_set=$2

  local concat_minus=concat_minus_$(basename $dissim_against_ref_fasta .fa).fa
  touch $concat_minus

  for fasta in "${dissim_set[@]}";
   do
    if [[ $fasta != $dissim_against_ref_fasta  ]]; then
      echo $fasta >> $concat_minus
    fi
   done
} # end createDissimilarConcatMinusFasta()

createFinalConcat(){
  local fa_to_concat=$1

  touch concat_transcriptome.fa

  for fasta in $fa_to_concat;
   do
     echo $fasta >> concat_transcriptome.fa
   done
} # end createFinalConcat

main $@
