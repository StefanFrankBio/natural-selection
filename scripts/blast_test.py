import argparse
import csv
import nstools


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--blast")
    parser.add_argument("-o", "--output")
    return parser.parse_args()


def main():
    args = parse_args()
    with open(args.blast, "r", ) as handle:
        csv_reader = csv.reader(handle, delimiter="\t")
        for blast_record in csv_reader:
            mismatch_count = int(blast_record[2])
            gap_count = int(blast_record[3])
            if min(mismatch_count, gap_count) > 0:
                aligned_reference = blast_record[5]
                aligned_variant = blast_record[4]
                absolute_pos = int(blast_record[1]) - 1
                variant, reference, insert_records = nstools.handle_inserts(aligned_reference, aligned_variant, absolute_pos)
                variant, deletion_records = nstools.handle_deletions(reference, variant, absolute_pos)
                variant, ambig_records = nstools.handle_ambiguities(reference, variant, absolute_pos)
                substitutions = nstools.find_substitutions(reference, variant, absolute_pos)
                variant_record = substitutions + insert_records + deletion_records + ambig_records
                variant_record = nstools.sort_by_element(variant_record, 0)
                nstools.write_seperated(f"{args.output}/{blast_record[0]}.tsv", variant_record)
                #nstools.vr_to_table("blast_test/variant_record.db", f"'{blast_record[0]}'", variant_record)


if __name__ == "__main__":
    main()
