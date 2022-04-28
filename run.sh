OUT_DIR=$1
THREADS=$2

TEST_LEN=1000
TEST_TYPE="d"
TEST_NUMBER=1000
TEST_VARIATIONS=10
REFERENCE="${OUT_DIR}/reference.fasta"
VARIANTS="${OUT_DIR}/variants.fasta"
TEST_VR_DIR="${OUT_DIR}/test_variant_records"

SPLIT_DIR="${OUT_DIR}/split_variants"
ALIGN_DIR="${OUT_DIR}/alignments"
VR_DIR="${OUT_DIR}/variant_records"

mkdir $TEST_VR_DIR
python scripts/test.py -l $TEST_LEN -t $TEST_TYPE -n $TEST_NUMBER -v $TEST_VARIATIONS -i -o $VARIANTS -r $REFERENCE -e $TEST_VR_DIR
seqkit split -i -2 -O $SPLIT_DIR $VARIANTS
rm "${VARIANTS}.seqkit.fai"
mkdir $ALIGN_DIR
find $SPLIT_DIR -type f | xargs -n 1 --max-procs=$THREADS stretcher -bsequence $REFERENCE -adirectory3 $ALIGN_DIR -aformat fasta -auto
rm -r $SPLIT_DIR
rm -r $VR_DIR
mkdir $VR_DIR
find $ALIGN_DIR -type f | xargs --max-args=1 --replace=1 --max-procs=$THREADS python scripts/variant_record.py --input 1 --reference $REFERENCE --output $VR_DIR
rm -r $ALIGN_DIR