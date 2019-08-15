
# blastn against ncbi_ise6 database, input query .fa as well as the output filepath in cmd line to script

# $1 : database
# $2 : query .fa
# $3 : output filename

blastn -outfmt 6 -num_threads 8 -db $1 -query $2 -out $3 
