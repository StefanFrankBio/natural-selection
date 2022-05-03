import argparse
import nstools

import sqlite3


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--reference")
    parser.add_argument("-v", "--variant")
    parser.add_argument("-d", "--database")
    parser.add_argument("-o", "--output")
    return parser.parse_args()


def main():
    args = parse_args()
    reference = nstools.read_fasta(args.reference)
    
    con = sqlite3.connect(args.database)
    con.row_factory = lambda cursor, row: row[0]
    cur = con.cursor()
    variants = cur.execute("SELECT name FROM sqlite_master where type='table'")
    for variant in variants:
        dNdS = [variant]
        for frame in range(3):
            codons = nstools.split_codons(str(reference.seq)[frame:])
            trans_table = nstools.build_trans_table()
            synonymous_sites = []
            for codon in codons:
                sns = nstools.single_substitutions(codon)
                translations = [trans_table[s] for s in sns]
                is_syn = [amino == trans_table[codon] for amino in translations]
                synonymous_sites += nstools.count_synonymous_sites(is_syn)
            variant_records = nstools.search_vr_table(args.database, variant, "Reference IN ('A','C','G','T') AND Variant IN ('A','C','G','T')")
            classification = nstools.classify_substitutions(variant_records, codons, frame, trans_table)
            dNdS.append(nstools.dNdS(classification, synonymous_sites))
        with open(args.output, "a") as handle:
            dNdS = "\t".join([str(i) for i in dNdS])
            print(dNdS, file=handle)


if __name__ == "__main__":
    main()
