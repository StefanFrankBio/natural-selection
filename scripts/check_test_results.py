import argparse

from numpy import var
import nstools


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--reference")
    parser.add_argument("-t", "--test")
    return parser.parse_args()

def main():
    args = parse_args()
    reference = nstools.read_fasta(args.reference)
    test_variant = nstools.read_fasta(args.test)
    variant_record = nstools.read_variant_records(f"../natural-selection-data/variant_records/{test_variant.id}")
    variant = nstools.reconstruct_variant(reference.seq, variant_record)
    if (variant == test_variant.seq) == False:
        print(f">{test_variant.id}")
        print(variant)

if __name__ == "__main__":
    main()