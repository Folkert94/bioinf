# Code Project group T20: Carlos, Folkert, Jurriaan, Lianne
# Edited code for Fundamentals of Bioinformatics Project

#!/usr/bin/python

import itertools
import argparse
import csv


def retrieve_scop_data(scop_file):
    """
    Reads a databaset file from SCOP and returns a dictionary with the protein IDS mapping.
    :param scop_file: database file containing mapping of PDB's to SCOP ID's.
    :return: dictionary containing the mapping of PDBs (keys) and their corresponding SCOP data (values).
    """

    scop_data = {}

    ##########################
    ### START CODING HERE ####
    ##########################
    # You can parse SCOP data in various ways. E.g. you can use dictionary of dictionaries
    # {proteinID: {"class": class, "fold": fold, "superfamily": superfamily, 'family': family}}
    scop_db = open(scop_file)
    for line in scop_db:
        if "#" in line.split()[0]:
            continue
        row = line.split()
        scop_data[row[1]] = [row[5]]

    # can possible contain just one pdb while there were more
    # print(scop_data["2bo9"])
    ########################
    ### END CODING HERE ####
    ########################
    return scop_data


def compute_similarity_score(prot1_scop, prot2_scop):
    """
    Computes the score for two proteins on the basis of the data from SCOP database.
    :param prot1_scop: data for protein 1 from SCOP database.
    :param prot2_scop: data for protein 2 from SCOP database.
    :return: similarity score (float)
    """
    ##########################
    ### START CODING HERE ####
    ##########################
    # You need to decide whether you need this function for SCOP database.
    pass


    ########################
    ### END CODING HERE ####
    ########################


def check_similarity_for_protein_pair(prot1_scop, prot2_scop):
    """
    Returns the similarity score between two proteins.
    :param prot1_scop: data for protein 1 from SCOP database.
    :param prot2_scop: data for protein 2 from SCOP database.
    :param pair: a tuple with the UniProt IDs of the two proteins to compare.
    :return: "different", "similar" or "ambiguous".
    """
    ##########################
    ### START CODING HERE ####
    ##########################
    data1 = prot1_scop.split(",")
    data2 = prot2_scop.split(",")
    prot1_f, prot1_sf, prot1_cf = data1[3].split("=")[1], data1[2].split("=")[1], data1[1].split("=")[1]
    prot2_f, prot2_sf, prot2_cf = data2[3].split("=")[1], data2[2].split("=")[1], data2[1].split("=")[1]
    if prot1_f == prot2_f:
        return "similar"
    if prot1_sf != prot2_sf and prot1_cf == prot2_cf:
        return "ambiguous"
    else:
        return "different"

    ########################
    ### END CODING HERE ####
    ########################

# If you will use the numeric score for SCOP (similar to GO), you may want to use check_similarity_for_protein_pair
# with other arguments. See the example below.
# def check_similarity_for_protein_pair(score, threshold):
#    pass

def generate_all_possible_protein_pairs(protein_ids):
    """
    Returns a list containing all unique protein pairs.
    :param protein_ids: list of all proteins IDs.
    :return: list of possible unique protein pairs.
    """
    pairs = list()
    ##########################
    ### START CODING HERE ####
    ##########################
    # You can add a pair of proteins to the list using the following code:
    # pairs.append((protein1, protein2))
    # Generate all possible combinations of IDs
    uniprot_id_list = []
    uniprot_ids = open(protein_ids)

    # Remove space in id_list
    for line in uniprot_ids:
        query = line.strip()
        uniprot_id_list.append(query)

    # Make possible combinations
    for pair in itertools.combinations(uniprot_id_list, 2):
        pairs.append(pair)

    # Remove duplicates
    pairs = set(pairs)
    pairs = set((a,b) if a<=b else (b,a) for a,b in pairs)

    ########################
    ### END CODING HERE ####
    ########################
    return list(pairs)


def assign_homology(scop_dict, protein_ids_pdbs, pairs):
    """
    Computes the similarity score between all protein pairs from the list, and decides if two proteins are homologs
    (different, ambiguous or similar).
    :param scop_dict: dictionary containing the mapping of PDBs (keys) and their corresponding SCOP data (values).
    :param protein_ids_pdbs: dictionary with UniprotID as key and PDB ID as a value.
    :param pairs: list of all possible unique protein pairs.
    :return: dictionary with UniProt ID (key), similarity(different, ambiguous or similar).
    """
    scop_homology = {}
    ##########################
    ### START CODING HERE ####
    ##########################
    # You should remember to take care about the proteins that are not in the SCOP database.

    for pair in pairs:
        pdb1 = protein_ids_pdbs[pair[0]][0].lower()
        pdb2 = protein_ids_pdbs[pair[1]][0].lower()
        try:
            sim_score = check_similarity_for_protein_pair(scop_dict[pdb1][0], scop_dict[pdb2][0])
            scop_homology[(pair[0], pair[1])] = sim_score
        except:
            continue
    ########################
    ### END CODING HERE ####
    ########################
    return scop_homology


def write_results(filename, scop_homology):
    """
    Writes in an output file the all of the protein pairs and their similarity/dissimilarity.
    :param output_filename: the name of the output file.
    :param scop_homology: dictionary (keys: protein pairs as tuples; values: one of the value - different/similar/ambiguous)

    """
    with open(filename, "w") as f:
        for (p1, p2), value in scop_homology.items():
            f.write("\t".join([p1, p2, value]) + "\n")


def read_protein_ids_file(filename):
    """
    Returns the list of UniProt IDs contained in the input file.
    :param filename:
    :return: list of UniProt IDs in the input file.
    """
    ##########################
    ### START CODING HERE ####
    ##########################
    uniprot_file = open(uniprot_filename, "r")
    uniprot_list = []
    for line in uniprot_file:
        uniprot_id = line.strip()
        # Fetch the fasta formatted sequence for each uniProtID.
        uniprot_list.append(uniprot_id)

    return uniprot_id_list
    #######################
    ### END CODING HERE ###
    #######################



def read_lookup_table(filename):
    """
    Reads the specified file and returns the dictionary with UniprotID as key and PDB ID as a value.
    :param filename: file with the mapping between Uniprot ids and PDB ids.
    :return: dictionary with UniprotID as key and PDB ID as a value.
    """
    ##########################
    ### START CODING HERE ####
    ##########################
    pdb_id_file = open(filename, "r")
    uniprot_pdb_dict = {}
    for line in pdb_id_file:
        pdb_id = str(line[7:-1])
        uniprot_id_for_dict = str(line[:-6])
        uniprot_pdb_dict.setdefault(uniprot_id_for_dict,[]).append(pdb_id)

    return uniprot_pdb_dict
    #######################
    ### END CODING HERE ###
    #######################



def main(uniprot_filename, output_filename, pdb_id_file, scop_file):
    ##########################
    ### START CODING HERE ####
    ##########################
    scop_dict = retrieve_scop_data(scop_file)
    protein_ids_pdbs = read_lookup_table(pdb_id_file)
    pairs = generate_all_possible_protein_pairs(uniprot_filename)
    scop_homo_dict = assign_homology(scop_dict, protein_ids_pdbs, pairs)
    write_results(output_filename, scop_homo_dict)
    #######################
    ### END CODING HERE ###
    #######################


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='The script retrieves data from SCOP database (from local file)'
                                                 ' and provides an output file'
                                                 ' with the strings "<id1>   <id2>   <similarity>", where'
                                                 ' <similarity> is a string with one of the values from '
                                                 ' different/ambiguous/similar.')

    parser.add_argument("-ids", "--uniprot_filename", help="File with the protein Uniprot IDs", required=True)
    parser.add_argument("-pdb", "--pdb_id_file", help="File with the mapping between Uniprot ids and PDB ids", required=True)
    parser.add_argument("-s", "--scop_file", help="SCOP database file", required=True)
    parser.add_argument("-o", "--output_filename", help="Output file name", required=True)

    args = parser.parse_args()

    uniprot_filename = args.uniprot_filename
    output_filename = args.output_filename
    pdb_id_file = args.pdb_id_file
    scop_file = args.scop_file

    main(uniprot_filename, output_filename, pdb_id_file, scop_file)
