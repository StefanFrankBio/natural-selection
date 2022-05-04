OUT_DIR=$1
THREADS=$2

REFERENCE="${OUT_DIR}/reference.fasta"
VARIANTS="${OUT_DIR}/variants.fasta"
BLAST_DIR="${OUT_DIR}/blast_results"
VR_DIR="${OUT_DIR}/variant_records"
VR_DB="${OUT_DIR}/variant.db"


rm -r $VR_DIR
mkdir $VR_DIR


makeblastdb -in $VARIANTS -dbtype nucl
blastn -query $REFERENCE -db $VARIANTS -num_threads=$THREADS -max_target_seqs=758775 -out blast_test/results.out -outfmt '6 sseqid qstart mismatch gapopen sseq qseq'
awk -F "\t" '{print $0 > "blast_test/blast_results/"($1)}' blast_test/results.out 
find $BLAST_DIR -type f | xargs --max-args=1 --replace=1 --max-procs=$THREADS python scripts/blast_test.py --blast 1 --output $VR_DIR
