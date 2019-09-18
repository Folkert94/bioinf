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
    i = 0
    for line in scop_db:
        if i > 3:
            row = line.split()
            temp = row[3].split(".")
            scop_data[row[1]] = {"class": temp[0], "fold": temp[1], "superfamily": temp[2], "family": temp[3]}
        i = i + 1

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
    pass


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

    # NOTE TO SELF
    # We ignore all the pairs that are not in the database as proposed by the assignment;
    # Did add them to the dictionary.


    for pair in pairs:
        # Match UniprotID to PDB ID
        pdb_ids_1 = protein_ids_pdbs[pair[0]]
        pdb_ids_2 = protein_ids_pdbs[pair[1]]
        # print(pdb_ids_1, pdb_ids_2)

        # Compute all pdb_id combinations
        pair_combinations = []
        for combination in itertools.product(pdb_ids_1, pdb_ids_2):
            pair_combinations.append(combination)

        # Check all against all in SCOP
        for pair_combination in pair_combinations:
            try:
                combi_1 = scop_dict[pair_combination[0].lower()]
            except:
                # print("pdb_id {0} not found".format(pair_combination[0]))
                scop_homology[(pair[0], pair[1])] = "not found"
                continue

            try:
                combi_2 = scop_dict[pair_combination[1].lower()]
            except:
                # print("pdb_id {0} not found".format(pair_combination[1]))
                scop_homology[(pair[0], pair[1])] = "not found"
                continue

            try:
                if combi_1["superfamily"] is combi_2["superfamily"] and combi_1["class"] is combi_2["class"] and combi_1["fold"] is combi_2["fold"]:
                    # print("{0} {1} similar".format(pair[0], pair[1]))
                    scop_homology[(pair[0], pair[1])] = "similar"
                if combi_1["fold"] is combi_2["fold"] and combi_1["superfamily"] is not combi_2["superfamily"]:
                    # print("{0} {1} ambiguous".format(pair[0], pair[1]))
                    scop_homology[(pair[0], pair[1])] = "ambiguous"
                else:
                    # print("{0} {1} different".format(pair[0], pair[1]))
                    scop_homology[(pair[0], pair[1])] = "different"
                # else:
                #     print("{0} AND {1} ARE UNREAL".format(pair_combination[0], pair_combination[1]))
                # print(combi_1["class"], combi_2["class"])
            except:
                pass

    ########################
    ### END CODING HERE ####
    ########################
    print(scop_homology)
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
