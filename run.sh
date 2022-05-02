OUT_DIR=$1
THREADS=$2

TEST_LEN=100
TEST_TYPE="d"
TEST_NUMBER=10
TEST_VARIATIONS=10
REFERENCE="${OUT_DIR}/reference.fasta"
VARIANT_DIR="${OUT_DIR}/test_variant_sequences"
TEST_VR_DIR="${OUT_DIR}/test_variant_records"
ALIGN_DIR="${OUT_DIR}/alignments"
VR_DIR="${OUT_DIR}/variant_records"

rm -r $VARIANT_DIR
rm -r $TEST_VR_DIR
rm -r $ALIGN_DIR
rm -r $VR_DIR
mkdir $VARIANT_DIR
mkdir $TEST_VR_DIR
mkdir $ALIGN_DIR
mkdir $VR_DIR

python scripts/generate_testdata.py -l $TEST_LEN -t $TEST_TYPE -n $TEST_NUMBER -v $TEST_VARIATIONS -i -o $VARIANT_DIR -r $REFERENCE -e $TEST_VR_DIR
find $VARIANT_DIR -type f | xargs -n 1 --max-procs=$THREADS stretcher -bsequence $REFERENCE -adirectory3 $ALIGN_DIR -aformat fasta -auto
find $ALIGN_DIR -type f | xargs --max-args=1 --replace=1 --max-procs=$THREADS python scripts/variant_record.py --input 1 --reference $REFERENCE --output $VR_DIR
find $VARIANT_DIR -type f | xargs --max-args=1 --replace=1 --max-procs=$THREADS python scripts/check_test_results.py --reference $REFERENCE --test 1

#KEEP FOR ONE-TIME EXECUTION
#seqkit split -i -2 -O $SPLIT_DIR $VARIANT_DIR
#rm "${VARIANT_DIR}.seqkit.fai"
#rm -r $SPLIT_DIR