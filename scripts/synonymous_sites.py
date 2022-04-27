import argparse
from Bio import SeqIO
import itertools
import pickle


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-c", "--codons")
    parser.add_argument("-s", "--sites")
    return parser.parse_args()


def read_fasta(filepath, multi=False):
    if multi == False:
        return next(SeqIO.parse(filepath, "fasta"))
    else:
        return list(SeqIO.parse(filepath, "fasta"))


def split_codons(sequence):
    overhang = len(sequence) % 3
    if overhang != 0:
        sequence = sequence[:-overhang]
    return [sequence[i:i+3] for i in range(0, len(sequence), 3)]


def single_substitutions(codon):
    substitutions = []
    for i, nt in enumerate(codon):
        possible_subs = ["A", "C", "G", "T"]
        possible_subs.remove(nt)
        for sub in possible_subs:
            substitutions.append(codon[:i] + sub + codon[i+1:])
    return substitutions


def build_trans_table():
    nucs = "TCAG"
    aminos = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    codons = (itertools.product(nucs, nucs, nucs))
    codons = ["".join(tpl) for tpl in codons]
    return dict(zip(codons, aminos))


def count_synonymous_sites(bool_list):
    return [sum(bool_list[x*3:(x+1)*3])/3 for x in range(3)]


def pickle_data(filepath, data):
    with open(filepath, "wb") as handle:
        pickle.dump(data, handle)


def main():
    args = parse_args()
    sequence_record = read_fasta(args.input)
    sequence = str(sequence_record.seq)
    codons = split_codons(sequence)
    trans_table = build_trans_table()
    synonymous_sites = []
    for codon in codons:
        sns = single_substitutions(codon)
        translations = [trans_table[s] for s in sns]
        is_syn = [amino == trans_table[codon] for amino in translations]
        synonymous_sites += count_synonymous_sites(is_syn)
    pickle_data(args.codons, codons)
    pickle_data(args.sites, synonymous_sites)


if __name__ == "__main__":
    main()