from Bio import SeqIO
import re
import itertools
import pickle
import random
from numpy import record

from pyparsing import line


def read_fasta(filepath: str, multi=False):
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


def count_synonymous_sites(bool_list):
    return [sum(bool_list[x*3:(x+1)*3])/3 for x in range(3)]


def handle_inserts(reference, sequence):
    sequence = list(sequence)
    insert_positions = [m.start() for m in re.finditer("-", reference)]
    insert_records = []
    for i in reversed(insert_positions):
        insert_records.append((i, "ins", sequence[i]))
        del sequence[i]
    return "".join(sequence), insert_records


def handle_deletions(reference, sequence):
    deletion_positions = [m.start() for m in re.finditer("-", sequence)]
    sequence = list(sequence)
    deletion_records = []
    for i in deletion_positions:
        deletion_records.append((i, reference[i], "del"))
        sequence[i] = reference[i]
    return "".join(sequence), deletion_records


def handle_ambiguities(reference, sequence):
    ambig_positions = list([m.start() for m in re.finditer("N", sequence)])
    sequence = list(sequence)
    ambig_records = []
    for i in ambig_positions:
        ambig_records.append((i, reference[i], "N"))
        sequence[i] = reference[i]
    return "".join(sequence), ambig_records


def pickle_data(filepath, data):
    with open(filepath, "wb") as handle:
        pickle.dump(data, handle)


def unpickle_data(filepath):
    with open(filepath, "rb") as handle:
        return pickle.load(handle)


def find_substitutions(sequence1, sequence2):
    xxx = zip(sequence1, sequence2)
    return [(i, n[0], n[1]) for i, n in enumerate(xxx) if n[0] != n[1]]


def build_trans_table():
    nucs = "TCAG"
    aminos = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    codons = (itertools.product(nucs, nucs, nucs))
    codons = ["".join(tpl) for tpl in codons]
    return dict(zip(codons, aminos))


def classify_substitutions(substitutions, codons, trans_table):
    classification = []
    for position, _, nucleotide in substitutions:
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


def write_seperated(filepath: str, args: list, seperator="\t") -> None:
    with open(filepath, "w") as handle:
        for row in args:
            print(*row, sep=seperator, file=handle)


def sort_by_element(data, idx):
    data.sort(key=lambda x: x[idx])
    return data


def generate_sequence(type: str, length: int) -> str:
    if type == "d":
        nucleotides = "ACGT"
        reference = "".join(random.choices(nucleotides, k=length))
    elif type == "p":
        aminoacids = "ARNDCQEGHILKMFPSTWYV"
        reference = "".join(random.choices(aminoacids, k=length))
    else:
        pass
    return reference


def test_variant_record(reference: str, type: str, length: int, variations: int, indels: bool) -> list:
    positions = sorted(random.sample(range(length), k=variations))
    ref_column = [reference[i] for i in positions]
    if type == "d":
        var_symbols = "ACGTN"
        var_column = random.choices(var_symbols, k=variations)
    elif type == "p":
        var_symbols = "ARNDCQEGHILKMFPSTWYVX"
        var_column = random.choices(var_symbols, k=variations)
    if indels == True:
        indel_count = random.randint(1, variations)
        indel_positions = random.sample(range(variations), k=indel_count)
        for i in indel_positions:
            ins_or_del = random.randint(0, 1)
            if ins_or_del == 0:
                ref_column[i] = "ins"
            else:
                var_column[i] = "del"
    else:
        pass
    variant_record = list(zip(positions, ref_column, var_column))
    return variant_record


def reconstruct_variant(reference: str, variant_record: list) -> str:
    variant = list(reference)
    for variation in variant_record[::-1]:
        position = variation[0]
        ref_symbol = variation[1]
        var_symbol = variation[2] 
        if ref_symbol == "ins":
            variant.insert(position, var_symbol)
        elif var_symbol == "del":
            del variant[position]
        else:
            variant[position] = var_symbol
    variant = "".join(variant)
    return variant


def read_variant_records(filepath: str) -> list:
    with open(filepath) as handle:
        variant_record = handle.readlines()
    variant_record = [line[:-1].split("\t") for line in variant_record]
    for variantion in variant_record:
        variantion[0] = int(variantion[0])
    return variant_record


if __name__ == "__main__":
    pass
