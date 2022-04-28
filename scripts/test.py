import argparse
import nstools
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--length", type=int)
    parser.add_argument("-t", "--type")
    parser.add_argument("-n", "--number", type=int)
    parser.add_argument("-v", "--variations", type=int)
    parser.add_argument("-i", "--indels", action=argparse.BooleanOptionalAction)
    parser.add_argument("-r", "--reference")
    parser.add_argument("-o", "--output")
    return parser.parse_args()

    
def main():
    args = parse_args()
    reference = nstools.test_reference(args.type, args.length)
    reference = SeqRecord(Seq(reference), id="test_reference", description="")
    with open(args.reference, "w") as handle:
        SeqIO.write(reference, handle, "fasta")
    with open(args.output, "w") as handle:
        for i in range(args.number):
            variant_record = nstools.test_variant_record(reference, args.type, args.length, args.variations, args.indels)
            variant = nstools.reconstruct_variant(reference, variant_record)
            variant = SeqRecord(Seq(variant), id=f"test_variant_{i}", description="")
            SeqIO.write(variant, handle, "fasta")


if __name__ == "__main__":
    main()