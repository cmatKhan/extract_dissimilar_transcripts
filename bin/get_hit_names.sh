
# get all unique values of first column of input
awk -v FS='\t' '{print $1}' $1 | uniq
