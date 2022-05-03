OUT_DIR=$1
THREADS=$2


REFERENCE="${OUT_DIR}/reference.fasta"
VARIANT_DIR="${OUT_DIR}/variant_sequences"
ALIGN_DIR="${OUT_DIR}/alignments"
VR_DIR="${OUT_DIR}/variant_records"
ERR_DIR="${OUT_DIR}/error_variant_records"


rm -r $ALIGN_DIR
rm -r $VR_DIR
rm -r $ERR_DIR
mkdir $ALIGN_DIR
mkdir $VR_DIR
mkdir $ERR_DIR


find $VARIANT_DIR -type f | xargs -n 1 --max-procs=$THREADS stretcher -bsequence $REFERENCE -adirectory3 $ALIGN_DIR -aformat fasta -auto
find $ALIGN_DIR -type f | xargs --max-args=1 --replace=1 --max-procs=$THREADS python scripts/variant_record.py --alignment 1 --reference $REFERENCE --variants $VARIANT_DIR --error $ERR_DIR --output $VR_DIR

ls $ERR_DIR | wc -l