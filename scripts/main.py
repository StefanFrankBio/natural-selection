import argparse
import nstools


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-c", "--codons")
    parser.add_argument("-s", "--sites")
    parser.add_argument("-r", "--reference")
    parser.add_argument("-o", "--output")
    return parser.parse_args()


def main():
    args = parse_args()
    reference = nstools.read_fasta(args.reference)
    alignment_records = nstools.read_fasta(args.input, True)
    aligned_reference = str(alignment_records[1].seq)
    aligned_sequence = str(alignment_records[0].seq)
    sequence, insert_records = nstools.handle_inserts(aligned_reference, aligned_sequence)
    sequence, deletion_records = nstools.handle_deletions(reference, sequence)
    sequence, ambig_records = nstools.handle_ambiguities(reference, sequence)
    codons = nstools.unpickle_data(args.codons)
    synonymous_sites = nstools.unpickle_data(args.sites)
    substitutions = nstools.find_substitutions(reference, sequence)
    trans_table = nstools.build_trans_table()
    synonymity = nstools.classify_substitutions(substitutions, codons, trans_table)
    nstools.write_seperated(
        f"{args.output}/{alignment_records[0].id}",
        alignment_records[0].id,
        substitutions,
        insert_records,
        deletion_records,
        ambig_records)


if __name__ == "__main__":
    main()