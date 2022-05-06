OUT_DIR=$1
THREADS=$2

REFERENCE="${OUT_DIR}/reference.fasta"
VARIANTS="${OUT_DIR}/variants.fasta"
BLAST_OUT="${OUT_DIR}/results.out"
BLAST_DIR="${OUT_DIR}/blast_results"
VR_DIR="${OUT_DIR}/variant_records"
VR_DB="${OUT_DIR}/variant.db"

rm $BLAST_OUT
rm -r $BLAST_DIR
rm -r $VR_DIR
rm $VR_DB
mkdir $BLAST_DIR
mkdir $VR_DIR

makeblastdb -in $VARIANTS -dbtype nucl
blastn -query $REFERENCE -db $VARIANTS -num_threads=$THREADS -max_target_seqs=758775 -out $BLAST_OUT -outfmt '6 sseqid qstart mismatch gapopen sseq qseq'
cat $BLAST_OUT | while read line
do
    echo "$line" | awk -v BLAST_DIR=$BLAST_DIR -F "\t" '{print $0 >> BLAST_DIR"/"($1)}'
done
find $BLAST_DIR -type f | xargs --max-args=1 --replace=1 --max-procs=$THREADS python scripts/blast_test.py --blast 1 --output $VR_DIR > test.txt

find $VR_DIR -type f | xargs --max-args=1 --replace=1 --max-procs=1 python scripts/vr_to_db.py --variant_record 1 --database $VR_DB >> test.txt