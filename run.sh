VARIANTS=$1
OUT_DIR=$2
REFERENCE=$3
THREADS=$4

SPLIT_DIR="${OUT_DIR}/split_variants"
ALIGN_DIR="${OUT_DIR}/alignments"
REF_CODONS="${OUT_DIR}/codons.pickle"
SYN_SITES="${OUT_DIR}/synonymous_sites.pickle"
VR_DIR="${OUT_DIR}/variant_records"

seqkit split -i -2 -O $SPLIT_DIR $VARIANTS
rm "${VARIANTS}.seqkit.fai"
mkdir $ALIGN_DIR
find $SPLIT_DIR -type f | xargs -n 1 --max-procs=$THREADS stretcher -bsequence $REFERENCE -adirectory3 $ALIGN_DIR -aformat fasta -auto
rm -r $SPLIT_DIR
python scripts/synonymous_sites.py --input $REFERENCE --codons $REF_CODONS --sites $SYN_SITES
rm -r $VR_DIR
mkdir $VR_DIR
find $ALIGN_DIR -type f | xargs --max-args=1 --replace=1 --max-procs=$THREADS python scripts/main.py --input 1 --reference $REFERENCE --output $VR_DIR
rm -r $ALIGN_DIR