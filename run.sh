VARIANTS=$1
OUT_DIR=$2
REFERENCE=$3
THREADS=$4

SPLIT_DIR="${OUT_DIR}/split_variants"
ALIGN_DIR="${OUT_DIR}/alignments"
REF_CODONS="${OUT_DIR}/codons.pickle"
SYN_SITES="${OUT_DIR}/synonymous_sites.pickle"
SUBS_DIR="${OUT_DIR}/substitutions"
SUBS_FILES="${SUBS_DIR}/*"

seqkit split -i -2 -O $SPLIT_DIR $VARIANTS
rm "${VARIANTS}.seqkit.fai"
mkdir $ALIGN_DIR
find $SPLIT_DIR -type f | xargs -n 1 --max-procs=$THREADS stretcher -bsequence $REFERENCE -adirectory3 $ALIGN_DIR -aformat fasta -auto
rm -r $SPLIT_DIR
python scripts/synonymous_sites.py --input $REFERENCE --codons $REF_CODONS --sites $SYN_SITES
mkdir $SUBS_DIR
find $ALIGN_DIR -type f | xargs --max-args=1 --replace=1 --max-procs=$THREADS python scripts/main.py --input 1 --reference $REFERENCE --codons $REF_CODONS --sites $SYN_SITES --output $SUBS_DIR
find $SUBS_DIR -type f -exec cat {} + > "${OUT_DIR}/all_substitutions.tsv"
rm -r $SUBS_DIR