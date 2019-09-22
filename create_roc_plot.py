#!/usr/bin/python

# This script reads and parses your previously obtained results out of a blast output file, and a benchmark output file.
# It then creates the corresponding ROC plot, exports the coordinates and calculates the AUC.

import argparse
import numpy
import matplotlib
matplotlib.use('AGG')
import pylab

def parse_blast_results(filename):
    """
    Parse every protein pair's e-value out of a BLAST results file.
    :param filename: input file with BLAST results.
    :return: dictionary with a tuple of two UniProt IDs (key), and the corresponding e-value from BLAST (value).
    """

    blast_results = {}

    with open(filename,'r') as f:
        for line in f:
            line = line.rstrip()
            arr = line.split("\t")

            if len(arr) != 3:
                print("Warning: the following line does not have three elements separated by a tab:\n", line)
            elif arr[0] == arr[1]:
                print("Warning: Comparing protein to itself:", arr[0])

            key = (arr[0], arr[1])
            if arr[2] == "NA":    # Substitute protein pairs whose e-value is
                value = 1e6       # not available with an e-value of 1 million
            else:
                value = float(arr[2])
            blast_results[key] = value

    return blast_results


def parse_benchmark_results(filename):
    """
    Parse every protein pair's classification out of the benchmark file.
    :param filename: input file with benchmark classifacations.
    :return: dictionary with a tuple of two UniProt IDs (key), and the corresponding call (value).
    """

    benchmark_results = {}

    with open(filename,'r') as f:
        for line in f:
            line = line.rstrip()
            arr = line.split("\t")

            if len(arr) < 3:
                print("Warning: the following line does not have three elements separated by a tab:\n", line)
            elif arr[0] == arr[1]:
                print("Warning: Comparing protein to itself:", arr[0])

            # Benchmark classifications are symmetric, so add both possible keys:
            key1 = (arr[0],arr[1])
            key2 = (arr[1],arr[0])
            value = arr[2]
            benchmark_results[key1] = value
            benchmark_results[key2] = value

    return benchmark_results


def integrate(x, y):
    """
    Calculate the Area Under the Curve (AUC) for a given list of coordinates
    :param x: a list of x-coordinates
    :param y: a list of y-coordinates
    :return: a float with the surface area under the curve described by x and y
    """

    auc = 0.
    last_x = x[0]
    last_y = y[0]
    for cur_x, cur_y in list(zip(x, y))[1:]:
        #########################
        ### START CODING HERE ###
        #########################
        auc += numpy.trapz([last_y, cur_y], [last_x, cur_x])

        #########################
        ###  END CODING HERE  ###
        #########################
        last_x = cur_x
        last_y = cur_y
    return auc


def roc_plot(blast_evalues, benchmark_dict, png_filename, evalue_blast, threshold):
    """
    Draw the ROC plot for a given set of e-values and corresponding benchmark classifications.

    :param blast_evalues: the dictionary produced by parse_blast_results()
    :param benchmark_dict: the dictionary produced by parse_benchmark_results()
    """

    ### Create the lists of coordinates

    x = [0] # array of the ROC plot's x-coordinates: False Positive Rate = FP/(FP+TN)
    y = [0] # array of the ROC plot's y-coordinates: True  Positive Rate = TP/(TP+FN)

    last_evalue = -1
    evalues = [(v, k) for k, v in blast_evalues.items()] # List of tuples consisting of (evalue, protein_pair)
    sorted_evalues = sorted(evalues)
    for evalue, protein_pair in sorted_evalues:

        #########################
        ### START CODING HERE ###
        #########################
        # Iterate through the protein pairs, in order of ascending e-value
        # Determine whether it is
        #    different -> actual negative, thus a false positive (x)
        #    similar   -> actual positive, thus a true positive (y)
        # Increase the respective value and add a new coordinate for every unique e-value
        # If the e-value is the same as the last one, only increase x or y of the last coordinate
        # Ignore entries in the benchmark_dict classified as "ambiguous" and decide how to handle blast NA results
        if evalue > last_evalue:
            x.append(x[-1])
            y.append(y[-1])
        if evalue <= threshold and benchmark_dict[protein_pair] == "different":
            x[-1] = x[-1] + 1
        if evalue <= threshold and benchmark_dict[protein_pair] == "similar":
            y[-1] = y[-1] + 1
            #########################
            ###  END CODING HERE  ###
            #########################
            last_evalue = evalue
    # In order to get the rates for every coordinate we divide by the total number (last entry)
    x = numpy.array(x) / float(x[-1])
    y = numpy.array(y) / float(y[-1])

    ### Figure out the AUC
    auc = integrate(x, y)

    ### Draw the plot and write it to a file
    pylab.plot(x, y, label="BLAST eval={0}, auc={1:.3f}".format(evalue_blast, auc))



    ### Write coordinates to a file
    with open(png_filename.split('.')[0] + '_xy.tsv','w') as f:
        for a,b in zip(x,y):
            f.write(str(a) + '\t' + str(b) + '\n')
    return auc


def roc_plot_psi(blast_evalues, benchmark_dict, png_filename, evalue_blast, threshold):
    """
    Draw the ROC plot for a given set of e-values and corresponding benchmark classifications.

    :param blast_evalues: the dictionary produced by parse_blast_results()
    :param benchmark_dict: the dictionary produced by parse_benchmark_results()
    """

    ### Create the lists of coordinates

    x = [0] # array of the ROC plot's x-coordinates: False Positive Rate = FP/(FP+TN)
    y = [0] # array of the ROC plot's y-coordinates: True  Positive Rate = TP/(TP+FN)

    last_evalue = -1
    evalues = [(v, k) for k, v in blast_evalues.items()] # List of tuples consisting of (evalue, protein_pair)
    sorted_evalues = sorted(evalues)
    for evalue, protein_pair in sorted_evalues:

        #########################
        ### START CODING HERE ###
        #########################
        # Iterate through the protein pairs, in order of ascending e-value
        # Determine whether it is
        #    different -> actual negative, thus a false positive (x)
        #    similar   -> actual positive, thus a true positive (y)
        # Increase the respective value and add a new coordinate for every unique e-value
        # If the e-value is the same as the last one, only increase x or y of the last coordinate
        # Ignore entries in the benchmark_dict classified as "ambiguous" and decide how to handle blast NA results
        if evalue > last_evalue:
            x.append(x[-1])
            y.append(y[-1])
        if evalue <= threshold and benchmark_dict[protein_pair] == "different":
            x[-1] = x[-1] + 1
        if evalue <= threshold and benchmark_dict[protein_pair] == "similar":
            y[-1] = y[-1] + 1
            #########################
            ###  END CODING HERE  ###
            #########################
            last_evalue = evalue
    # In order to get the rates for every coordinate we divide by the total number (last entry)
    x = numpy.array(x) / float(x[-1])
    y = numpy.array(y) / float(y[-1])

    ### Figure out the AUC
    auc = integrate(x, y)

    ### Draw the plot and write it to a file
    pylab.plot(x, y, label="PSIBLAST eval={0}, auc={1:.3f}".format(evalue_blast, auc), linestyle="dashed")



    ### Write coordinates to a file
    with open(png_filename.split('.')[0] + '_xy.tsv','w') as f:
        for a,b in zip(x,y):
            f.write(str(a) + '\t' + str(b) + '\n')
    return auc


def main(blast_results_map, benchmark_results_file, png_file):
    pylab.figure()
    pylab.plot([0,1],[0,1],'--k')
    pylab.xlabel('False Positive Rate')
    pylab.ylabel('True Positive Rate')
    pylab.title('Plots BLAST&PSI-Blast')


    evalues = [1.0, 10.0, 100.0]
    auc_list = []
    for i in range(len(evalues)):
        evalue_blast = evalues[i]
        # Parse the input files and retrieve every protein pair's e-value and benchmark classification.
        blast_evalues = parse_blast_results("{0}/{0}{1}".format(blast_results_map, i + 3))
        benchmark_results = parse_benchmark_results(benchmark_results_file)

        # Draw and save the ROC plot
        roc_plot(blast_evalues, benchmark_results, png_file, evalue_blast, evalues[i])

    for i in range(len(evalues)):
        evalue_psiblast = evalues[i]
        # Parse the input files and retrieve every protein pair's e-value and benchmark classification.
        blast_evalues = parse_blast_results("output_psiblast/output_psiblast{0}".format(i + 3))
        benchmark_results = parse_benchmark_results(benchmark_results_file)

        # Draw and save the ROC plot
        roc_plot_psi(blast_evalues, benchmark_results, png_file, evalue_psiblast, evalues[i])





    pylab.legend()
    pylab.savefig(png_file)



    # ADDITIONAL CODE FOR AUC VS EVALUE PLOTS
    # pylab.close()
    #
    # pylab.figure()
    # pylab.xlabel('log E-values')
    # pylab.ylabel('AUC')
    # pylab.title('AUC vs E-values PSI-BLAST')
    # pylab.plot(numpy.log(evalues), auc_list)
    # pylab.savefig("auc-vs-evalue PSI-BLAST.png")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Draw and save a ROC plot to a file")
    parser.add_argument("-iblast","--input_blast_results", help="map with tab-separated (PSI-)BLAST results files", required=True)
    parser.add_argument("-ibench","--input_benchmark_results", help="tab-separated benchmark classification file", required=True)
    parser.add_argument("-o", "--output_png", help="output png file", required=True)

    args = parser.parse_args()

    blast_results_map = args.input_blast_results
    benchmark_file = args.input_benchmark_results
    png_file = args.output_png

    main(blast_results_map,benchmark_file, png_file)
