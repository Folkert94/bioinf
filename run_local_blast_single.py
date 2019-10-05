
# Code Project group T20: Carlos, Folkert, Jurriaan, Lianne
# Edited code for Fundamentals of Bioinformatics Project
# This script runs the local version of (PSI-)BLAST.

#!/usr/bin/python
import argparse
import subprocess
import math
import numpy as np
import itertools
import matplotlib
matplotlib.use('AGG')
import pylab


def blast(db, query, evalue, psiblast=False, query_folder="./queries/"):
    """
    This function executes blast or psi-blast for the given query and db.
    :param db: database filename
    :param query: query filename
    :param query_folder: query folder name
    :param psiblast: True if PSI-BLAST should be used; False for normal BLAST
    :return: result from blast run
    """
    if psiblast:
        cmd = "psiblast -query {0}{1}.fasta -db {2} -outfmt '6 qacc sacc evalue' -num_iterations 3 -evalue {3}".format(query_folder, query, db, evalue)

    else:
        cmd = "blastp -query {0}{1}.fasta -db {2} -outfmt '6 qacc sacc evalue' -evalue {3}".format(query_folder, query, db, evalue)

    # Running shell command in python script. See https://docs.python.org/2/library/subprocess.html#popen-constructor
    p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         close_fds=True)
    blast_result = p.stdout.read().decode("utf8")
    return blast_result


def parse_blast_result(blast_result, blast_dict):
    """
    This function parses the output of (PSI-)BLAST and stores the result in blast_dict (defined in main()).
    :param blast_result: output  obtained after running (PSI-)BLAST
    :param blast_dict: dictionary where the sequence's alignments will be stored [Keys: (id1,id2) Values: E-values]
    :return: dictionary storing protein pair as tuple (keys) and the corresponding e-value (values)
    """
    for line in blast_result.split("\n"):
        if line and line[0] != "#" and line[0] != "[":
            try:
                splitted_line = line.split()
                query = splitted_line[0].split("|")[1]
                subject = splitted_line[1].split("|")[1]

                pair = (query, subject)
                result = float(splitted_line[2])
                blast_dict[pair] = result

            except IndexError:
                if not line.endswith('CONVERGED!'):
                    print ("\tCould not parse (psi-)blast response line:\n\t"+ line)

    return blast_dict


def write_output(uniprot_ids, output_filename, blast_dict):
    """
    This function writes the scores of all-against-all protein pairs to the output file.
    :param blast_result: output result obtained after running (PSI-)BLAST.
    :param blast_dict: dictionary where the sequence's alignments will be stored [Keys: (id1,id2) Values: E-values].
    :return: dictionary storing protein pair as tuple (keys) and the corresponding e-value (values).
    """
    with open(output_filename, "w") as f:
        # Generate all possible combinations of IDs
        combinations_ids = itertools.product(uniprot_ids, uniprot_ids)
        for pair in combinations_ids:
            # Check if the pair has unique elements
            if len(pair) == len(set(pair)):
                pair_str = "\t".join(pair)
                if pair in blast_dict:
                    f.write(pair_str + "\t" + str(blast_dict[pair]) + "\n")
                else:
                    f.write(pair_str + "\t" + "NA\n")


def plot_evalue_distribution(blast_dict, png_filename, evalue=100000):
    """
    This function plots the distribution of the log(e-value). pseudocount is added to avoid log(0).
    The pseudocount in this case is the smallest non-zero e-value divided by 1000.
    :param blast_dict: dictionary containing all the sequences alignments and e-values [Keys: (id1,id2) Values: E-values]
    :param png_filename: png output file to save the distribution plot.
    :param evalue: threshold for e-value. If no threshold specified, arbitrary 100000 will be used.

    """
    sorted_e_val = sorted(blast_dict.values())
    nonzero_indices = np.nonzero(sorted_e_val)[0]
    pseudo_count = sorted_e_val[nonzero_indices[0]] / 1000.0
    pylab.hist(list(map(lambda x: math.log10(x + pseudo_count), blast_dict.values())))
    pylab.xlabel("log(e-value)")
    pylab.ylabel("Frequency")
    pylab.savefig(png_filename)

    total_zeros = sum(1 for i in blast_dict.values() if i < evalue)
    print("The number of e-values lower than the threshold ({0}) is: {1}".format(evalue, total_zeros))


def main(uniprot_id_file, query_folder, db, evalue, psiblast, output_filename, output_png):

    # The blast_dict dictionary will be used to store protein pair and the corresponding e-value.
    # Keys for blast_dict are the combination of query and subject/hit, e.g.:
    # key             = (query, subject)
    # blast_dict[key] = e_value
    blast_dict = {}
    # uniprot_id_list is a list to store all UniProt IDs contained in uniprot_id_file.
    uniprot_id_list = []
    print("Running blast for e-value {0}.".format(evalue))
    uniprot_ids = open(uniprot_id_file)
    for line in uniprot_ids:
        query = line.strip()

        blast_result = blast(db, query, evalue, psiblast)
        uniprot_id_list.append(query)
        parse_blast_result(blast_result, blast_dict)

    output_filename1 = "{0}{1}".format(output_filename, evalue)
    uniprot_ids.close()
    write_output(uniprot_id_list, output_filename1, blast_dict)
    plot_evalue_distribution(blast_dict, "{0}{1}.png".format(output_png, evalue), evalue)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automatically running BLAST and PSI-BLAST")
    parser.add_argument("-ids", "--protein_ids_file", help="the list of UniProt IDs", required=True)
    parser.add_argument("-q", "--query_folder", help="the query folder", required=True)
    parser.add_argument("-db", "--db", help="the fasta file of the database", required=True)
    parser.add_argument("-eval", "--evalue", help="the desired e-value, default is 10", default=10, required=False)
    parser.add_argument("-o", "--output_file", help="output file", required=True)
    parser.add_argument("-opng", "--output_png", help="output png file", default="DistributionEValue.png", required=False)
    parser.add_argument("-psi", "--psiblast", dest="psiblast", action="store_true", help="If flagged, run PSI-BLAST instead of BLASTP")

    args = parser.parse_args()

    # Assign the parsed arguments to the corresponding variables.
    uniprot_id_file = args.protein_ids_file
    query_folder = args.query_folder
    db = args.db
    evalue =args.evalue
    psiblast = args.psiblast # True or False
    output_filename = args.output_file
    output_png = args.output_png

    main(uniprot_id_file, query_folder, db, evalue, psiblast, output_filename, output_png)
