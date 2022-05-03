import argparse
import nstools


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--alignment")
    parser.add_argument("-r", "--reference")
    parser.add_argument("-v", "--variants")
    parser.add_argument("-o", "--output")
    parser.add_argument("-e", "--error")
    return parser.parse_args()


def main():
    args = parse_args()
    reference = nstools.read_fasta(args.reference)
    alignment_records = nstools.read_fasta(args.alignment, True)
    aligned_reference = str(alignment_records[1].seq)
    aligned_variant = str(alignment_records[0].seq)
    variant, insert_records = nstools.handle_inserts(aligned_reference, aligned_variant)
    variant, deletion_records = nstools.handle_deletions(reference, variant)
    variant, ambig_records = nstools.handle_ambiguities(reference, variant)
    substitutions = nstools.find_substitutions(reference, variant)
    variant_record = substitutions + insert_records + deletion_records + ambig_records
    variant_record = nstools.sort_by_element(variant_record, 0)
    reconstructed_variant = nstools.reconstruct_variant(reference.seq, variant_record)
    variant = nstools.read_fasta(f"{args.variants}/{alignment_records[0].id}.fasta")
    if variant.seq == reconstructed_variant:
        nstools.write_seperated(f"{args.output}/{alignment_records[0].id}", variant_record)
    else:
        nstools.write_seperated(f"{args.error}/{alignment_records[0].id}", variant_record)


if __name__ == "__main__":
    main()