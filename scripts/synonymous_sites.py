import argparse
import nstools


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-c", "--codons")
    parser.add_argument("-s", "--sites")
    return parser.parse_args()


def main():
    args = parse_args()
    sequence_record = nstools.read_fasta(args.input)
    sequence = str(sequence_record.seq)
    codons = nstools.split_codons(sequence)
    trans_table = nstools.build_trans_table()
    synonymous_sites = []
    for codon in codons:
        sns = nstools.single_substitutions(codon)
        translations = [trans_table[s] for s in sns]
        is_syn = [amino == trans_table[codon] for amino in translations]
        synonymous_sites += nstools.count_synonymous_sites(is_syn)
    nstools.pickle_data(args.codons, codons)
    nstools.pickle_data(args.sites, synonymous_sites)


if __name__ == "__main__":
    main()