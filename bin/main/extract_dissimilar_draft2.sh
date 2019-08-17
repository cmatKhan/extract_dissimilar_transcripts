#########################################################################################################################
# name: extract_dissimilar_main
# purpose: main script driving extract_dissimilar_transcripts
# input: Run script from location at which you wish to deposit the results directory
#        $1 : reference transcriptome
#        $2 : a directory containing (and only containing) the transcriptomes to be added (excluding reference transcriptome).
#             All transcriptomes must have extension .fa (easily taken care of if you do use symbollic links as you can rename the symbollic link rather than the raw data file)
#             If you store these directories separately on your local machine, create a
#              directory with symbollic links to the appropriate transcriptomes
#        $3 : project name
#        $4 : path to repo/bin
# output: a concatenated transcriptome consisting of the reference transcriptome +
#         the most dissimilar sequences from each of the other transcriptomes
# written by: chase mateusiak, chase.mateusiak@gmail.com, chase.mateusiak@ucsf.edu
# date: 20190814, edit 20190815
#
# Dependencies: blast+ > ___v, pandas etc from .py scripts(s), python3;
#               PLEASE NOTE: this is intended to be run from wynton qb3-dev3 (or any of the other dev nodes, but this one has the most resources)
#                            This will also work on your local machine if you have blast+ installed and blastn and makeblastdb are in your path.
#                            However, you must remove the module load lines from the makeBlastDB and makeDissimilarFasta methods
#
# Note: .fa must be standard -- pacbio names transcriptome needed to be manipulated to
#########################################################################################################################

# main ~/data/iscap_transcriptomes/raw_transcriptomes/ise6_rna_ncbi_20190808.fa ~/data/iscap_transcriptomes/add_transcriptomes ncbi_ise6_plus_all ~/tick/extract_dissimilar_transcripts/bin

main(){
  # rename cmd line input for clarity
  local reference_transcriptome=$1
  local input_transcriptomes=$2
  local project_name=$3
  local repo_bin=$4

  # create project directory at current location, cd, store path
  mkdirCd extract_dissimilar_${project_name}
  local project_dir=$(pwd)
  mkdir logs
  local logs_dir=${project_dir}/logs

  # make database from reference_transcriptome
  mkdirCd reference_database
  printf "\ncreate blast database from reference transcriptome"
  makeBlastDB $reference_transcriptome
  # store path to reference db (/path/to/dir/basename_of_database_files.*)
  local reference_db_path=$(realpath $(ls .))
  local reference_transcriptome_db=${reference_db_path}/$(basename $reference_db_path)

  # return to project directory
  cd $project_dir

  # create dir input_transcriptomes which will hold data for each input_transcriptome, cd, store path
  mkdirCd input_transcriptomes
  local input_transcriptomes_dir=$(pwd)

  # loop through input_transcriptomes, blast transcript against ref, get unique names of transcripts found to be similar to reference,
  # use unique names to extract transcripts without any similarity to reference and create .fa of these
  for transcriptome in $input_transcriptomes/*;
    do
      printf "\nbegin processing $transcriptome against reference transcriptome database\n"
      # create directory and subdirectories for each input_transcriptome with subdir blast_against_ref, cd into it
      mkdirCd $(basename $transcriptome .fa)
      mkdirCd blast_against_reference
      # extract transcripts without similarities to ref_transcriptome, output as .fa
      makeDissimilarFasta $reference_transcriptome_db $transcriptome $repo_bin
      # add .fa to dissimilar_to_ref_set
      dissimilar_fa=$(realpath $(find . -name "*.fa"))
      echo $dissimilar_fa >> ${project_dir}/logs/dissimilar_to_ref_set.txt
      # move out of individual input_transcriptome to input_transcriptome dir
      cd $input_transcriptomes_dir
    done

    # a list to store the .fasta of subset of transcripts from input_transcriptomes dissimilar from reference_transcriptome
    dissimilar_to_ref_set=$(createList ${project_dir}/logs/dissimilar_to_ref_set.txt)

    printf "\n\n\n\nthe set of fastas dissimilar to the reference transcriptome are ${dissimilar_to_ref_set}\n\n\n\n"
    printf "\n"

  # loop through each new_transcriptome dir, create concat_transcriptome minus the
  # individual new transcriptome, create new blast database from this concatMinus db
  # blastn the dissimilar .fa against concatMinus_db and then extract those dissimilar transcriptomes,
  # which will be the unique transcripts from the given new_transcriptome to be added to the ref_transcriptome
  # in the final concatenation

  for input_transcriptome_dir in $input_transcriptomes_dir/*;
    do
      printf "\nbegin processing $input_transcriptome_dir to extract .fa of transcripts dissimilar to other transcripts dissimilar to reference transcriptome\n"

      cd $input_transcriptome_dir
      # make new directory to hold data related to within_dissimilar set comparison
      mkdirCd blast_against_dissimilar

      # current input_transcriptome dissimilar_to_reference .fa
      local current_fasta=$(realpath $(find ../blast_against_reference -name "*.fa"))

      # input dissimilar .fa in input_transcriptome_dir and array of all dissimilar .fa, create concat .fa of all dissimilar_input_transcriptomes EXCEPT current input_transcriptome
      mkdirCd concat_minus_fa

      printf "\ninput to createDissimilar, remove_fasta: $current_fasta\n set:$dissimilar_to_ref_set\n"

      # create concat transcriptome minus current_fasta and a filename for new .fa
      local current_fasta_bn=$(basename $1 .fa)
      local concat_minus=concat_minus_${current_fasta_bn}.fa

      for fa in $dissimilar_to_ref_set;
       do
        printf "\ncreateDissimilar with fasta from set: $fa, current_fasta: $current_fasta\n"
        if [[ $fa != $current_fasta ]]; then
          printf "\nconcat $fa to $concat_minus\n"
          cat $fa >> $concat_minus
        fi
       done
      cd ..
      # make blast database from concatMinus .fa
      printf "\ninput to makeBlastDB $(realpath $(find ./concat_minus_fa -name "*.fa"))\n"
      makeBlastDB $(realpath $(find ./concat_minus_fa -name "*.fa"))
      # store path to concatMinus database
      local databaseMinus=$(realpath $(find . -name "*_db"))
      local databaseMinus_name=$(basename $databaseMinus)
      local databaseMinus_path=${databaseMinus}/${databaseMinus_name}
      # create .fa of transcripts which dissimilar to the transcripts which are also dissimilar to the reference transcriptome
      printf "\npath to database is $databaseMinus_path\n"
      makeDissimilarFasta $databaseMinus_path $dissimilar_against_reference_fasta $repo_bin
      # add this .fa to array which stores final .fa to concat to reference_transcriptome
      final_concat_fa=$(realpath $(find . -name "*.fa"))
      echo $final_concat_fa >> ${project_dir}/logs/dissimilar_dissimilar_fa_to_concat.txt
      cd $input_transcriptomes_dir
    done

    transcripts_to_concat=$(createList ${project_dir}/logs/dissimilar_dissimilar_fa_to_concat.txt)

  cd $project_dir
  createFinalConcat $ref_transcriptome $transcripts_to_concat
  exit 0
} # end main()

mkdirCd(){
  mkdir $1
  cd $1
} # end mkdirCd

makeBlastDB (){
  module load Sali
  module load blast+

  local fasta=$1
  local name=$(basename $fasta .fa)
  local db_name=${name}_db
  mkdirCd ${db_name}
  makeblastdb -parse_seqids -dbtype nucl -in $fasta -out $db_name
  cd ..
} # end makeBlastDB

makeDissimilarFasta(){
  module load Sali
  module load blast+
  module load anaconda

  local CONDA_BASE=$(conda info --base)
  source ${CONDA_BASE}/etc/profile.d/conda.sh
  # activate tae1 in atanas' wynton acct b/c it is the Bio modules
  conda activate tae1

  local db=$1
  local query_fasta=$2
  local query_fasta_bn=$(basename $2 .fa)
  local output_tsv=${query_fasta_bn}.tsv
  local output_noext=$(basename $output_tsv .tsv)
  local output_pairwise=${output_noext}.pairwise
  local bin=$3

  printf "\nblast $query_fasta_bn , create .tsv and .pairwise\n"
  # blast query_fasta against db -- output both .tsv and pairwise comparisons
  blastn -outfmt 6 -num_threads 8 -db $db -query $query_fasta -out $output_tsv
  blastn -outfmt 0 -num_threads 8 -db $db -query $query_fasta -out $output_pairwise

  printf "\nextract unique transcript names of sequences with similarities to db from $output_noext\n"
  #extract unique names from output of above
  awk -v FS='\t' '{print $1}' $output_tsv | uniq > ${output_noext}_unique.tsv

  printf "\ncreate $2 of sequences without significant matches in db\n"
  # create .fa of transcripts from $query_fasta not found to have any similarity to transcripts in $db
  if [[ $query_fasta_bn == p* ]]; then
    #python3 /wynton/home/choulab/atanasdradkov/tick/extract_dissimilar_transcripts/bin/create_fa_blast_unmatched_pacbio.py $query_fasta ${output_noext}_unique.tsv
    python3 ${bin}/create_fa_blast_unmatched_pacbio.py $query_fasta ${output_noext}_unique.tsv
  else
    #python3 /wynton/home/choulab/atanasdradkov/tick/extract_dissimilar_transcripts/bin/create_fa_blast_unmatched.py $query_fasta ${output_noext}_unique.tsv
    python3 ${bin}/create_fa_blast_unmatched.py $query_fasta ${output_noext}_unique.tsv
  fi
  conda deactivate
} # end makeDissimilarFa()

# createDissimilarConcatMinusFasta(){
#   local current_fasta=$1
#   local current_fasta_bn=$(basename $1 .fa)
#   local dissim_set=$(echo $2)
#
#   local concat_minus=concat_minus_${current_fasta_bn}.fa
#   for fa in $dissim_set;
#    do
#     printf "\ncreateDissimilar with fasta from set: $fasta, current_fasta: $current_fasta\n"
#     if [[ $fa != $current_fasta ]]; then
#       printf "\nconcat $fasta to $concat_minus\n"
#       cat $fa >> $concat_minus
#     fi
#    done
# } # end createDissimilarConcatMinusFasta()

createFinalConcat(){
  printf "\ncreate final concat transcriptome\n"

  local fa_to_concat=$1
  local ref=$2

  touch dissimilar_transcripts.fa
  touch concat_ref_dissimilar.fa

  for fasta in $fa_to_concat;
   do
     printf "is this just going forever? why isn't it stopping?"
     cat $fasta >> dissimiliar_transcripts.fa
   done
   cat $2 >> concat_ref_dissimilar.fa
   cat dissimiliar_transcripts.fa >> concat_ref_dissimilar.fa
} # end createFinalConcat

createList(){
  file_to_arr=$1
  declare -a arr=()
  local length=$(wc -l $file_to_arr | cut -f1 -d ' ')
  for i in $(seq 1 $length);
   do
     awk "NR==$i {print}" $file_to_arr
   done
} # end createList()

main $@
