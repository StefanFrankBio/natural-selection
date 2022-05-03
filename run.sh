OUT_DIR=$1
THREADS=$2


REFERENCE="${OUT_DIR}/reference.fasta"
VARIANT_DIR="${OUT_DIR}/variant_sequences"
ALIGN_DIR="${OUT_DIR}/alignments"
ERR_DIR="${OUT_DIR}/error_variant_records"
VR_DB="${OUT_DIR}/variant.db"
DNDS="${OUT_DIR}/dNdS.tsv"


rm $VR_DB
rm $DNDS
rm -r $ALIGN_DIR
rm -r $ERR_DIR
mkdir $ALIGN_DIR
mkdir $ERR_DIR


find $VARIANT_DIR -type f | xargs -n 1 --max-procs=$THREADS stretcher -bsequence $REFERENCE -adirectory3 $ALIGN_DIR -aformat fasta -auto -verbose
find $ALIGN_DIR -type f | xargs --max-args=1 --replace=1 --max-procs=$THREADS python scripts/variant_record.py --alignment 1 --reference $REFERENCE --variants $VARIANT_DIR --database $VR_DB --error $ERR_DIR
python scripts/dNdS.py -r $REFERENCE -v 1 -d $VR_DB -o $DNDS

ls $ERR_DIR | wc -l