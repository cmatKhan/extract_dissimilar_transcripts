
# make blast database from inputted .fasta

# $1 : input .fasta

# $2 : path to output file

makeblastdb -parse_seqids -dbtype nucl -in $1 -out $2

