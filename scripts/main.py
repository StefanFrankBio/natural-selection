import argparse
from Bio import SeqIO
from re import finditer, sub
import pickle
import itertools


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-c", "--codons")
    parser.add_argument("-s", "--sites")
    parser.add_argument("-r", "--reference")
    parser.add_argument("-o", "--output")
    return parser.parse_args()


def read_fasta(filepath, multi=False):
    if multi == False:
        return next(SeqIO.parse(filepath, "fasta"))
    else:
        return list(SeqIO.parse(filepath, "fasta"))


def handle_inserts(reference, sequence):
    sequence = list(sequence)
    insert_positions = [m.start() for m in finditer("-", reference)]
    insert_records = []
    for i in reversed(insert_positions):
        insert_records.append((i,sequence[i]))
        del sequence[i]
    return "".join(sequence), insert_records


def handle_deletions(reference, sequence):
    deletion_positions = [m.start() for m in finditer("-", sequence)]
    sequence = list(sequence)
    deletion_records = []
    for i in deletion_positions:
        deletion_records.append((reference[i],i))
        sequence[i] = reference[i]
    return "".join(sequence), deletion_records


def handle_ambiguities(reference, sequence):
    ambig_records = list([m.start() for m in finditer("N", sequence)])
    sequence = list(sequence)
    for i in ambig_records:
        sequence[i] = reference[i]
    return "".join(sequence), ambig_records


def unpickle_data(filepath):
    with open(filepath, "rb") as handle:
        return pickle.load(handle)


def find_substitutions(sequence1, sequence2):
    xxx = zip(sequence1, sequence2)
    return [(n[0], i, n[1]) for i, n in enumerate(xxx) if n[0] != n[1]]


def build_trans_table():
    nucs = "TCAG"
    aminos = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    codons = (itertools.product(nucs, nucs, nucs))
    codons = ["".join(tpl) for tpl in codons]
    return dict(zip(codons, aminos))


def classify_substitutions(substitutions, codons, trans_table):
    classification = []
    for _, position, nucleotide in substitutions:
        codon_idx = position // 3
        if codon_idx < len(codons):
            original_codon = codons[codon_idx]
            position_in_codon = position % 3
            substituted_codon = original_codon[:position_in_codon] + nucleotide + original_codon[position_in_codon+1:]
            classification.append(trans_table[original_codon] == trans_table[substituted_codon])
    return classification


def dNdS(synonymity, synonymous_sites):
    syn_site_count = sum(synonymous_sites)
    nonsyn_site_count = len(synonymous_sites) - syn_site_count
    syn_subs = sum(synonymity)
    nonsyn_subs = len(synonymity) - syn_subs
    if syn_subs == 0:
        return "div by 0"
    else:
        return (nonsyn_subs/nonsyn_site_count)/(syn_subs/syn_site_count)


def sqlite_test(id, dNdS_ratio):
    import sqlite3
    con = sqlite3.connect("db_test.db")
    cur = con.cursor()
    cur.execute(f"UPDATE metadata SET dNdS = '{dNdS_ratio}' WHERE ID = '{id}'")
    con.commit()
    con.close()


def main():
    args = parse_args()
    reference = read_fasta(args.reference)
    alignment_records = read_fasta(args.input, True)
    aligned_reference = str(alignment_records[1].seq)
    aligned_sequence = str(alignment_records[0].seq)
    sequence, insert_records = handle_inserts(aligned_reference, aligned_sequence)
    sequence, deletion_records = handle_deletions(reference, sequence)
    sequence, ambig_records = handle_ambiguities(reference, sequence)
    codons = unpickle_data(args.codons)
    synonymous_sites = unpickle_data(args.sites)
    substitutions = find_substitutions(reference, sequence)
    trans_table = build_trans_table()
    synonymity = classify_substitutions(substitutions, codons, trans_table)
    dNdS_ratio = dNdS(synonymity, synonymous_sites)

    with open(f"{args.output}/{alignment_records[0].id}", "w") as handle:
        print(alignment_records[0].id, substitutions, insert_records, deletion_records, ambig_records, sep="\t", file=handle)

if __name__ == "__main__":
    main()