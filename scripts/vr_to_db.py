import argparse
import csv

from numpy import var
import nstools

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--variant_record")
    parser.add_argument("-d", "--database")
    return parser.parse_args()


def main():
    args = parse_args()
    seq_id = args.variant_record.split("/")[-1][:-4]
    with open(f"{args.variant_record}") as handle:
        csv_reader = csv.reader(handle, delimiter="\t")
        variant_record = [record for record in csv_reader]

    nstools.vr_to_table(args.database, f"'{seq_id}'", variant_record)
    print("A")


if __name__ == "__main__":
    main()
