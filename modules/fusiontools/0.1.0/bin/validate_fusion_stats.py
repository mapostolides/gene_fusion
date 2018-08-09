#! /usr/bin/env python

import sys
import argparse

# instantiate parser object
parser = argparse.ArgumentParser()

parser.add_argument('validated_fusions_file', action='store', help='A file containing a set of validated fusions', type=file)
parser.add_argument('fusion_cluster_file', action='store', help='Fusion reannation file clustered by head/tail genes. (.reann.cluster file)')

args = parser.parse_args()

def compare_validated_and_output_fusions(validated_fusions, output_fusions):
    """
    Compares the fusion pairs in the validated fusion file with the detected fusions in the output file
    Outputs results/statistics to an output file
    :param validated_fusions: a text file containing fusion gene pairs we would expect to be detected 
                              by the pipeline with the given input sequence data
    :param output_fusions: file merged.cff.reann.dnasupp.bwafilter.30.cluster, outputted by the pipeline
    """
    for line_val in validated_fusions:
        gene1_val, gene2_val = line_val.split()[1], line_val.split()[2]
        #print "Validated genes: {gene1_val} {gene2_val} ".format(gene1_val=gene1_val, gene2_val=gene2_val)

        # need to open this on every iteration, since once a file is open, it can only be iterated through once, unless opened again
        for line_out in open(output_fusions, 'r'):
            gene1_out, gene2_out = line_out.split()[1], line_out.split()[2]
            #print "output genes: {gene1_out} {gene2_out} ".format(gene1_out=gene1_out, gene2_out=gene2_out) 
            if (gene1_val == gene1_out) and (gene2_val == gene2_out):
                print gene1_val, gene2_val


compare_validated_and_output_fusions(args.validated_fusions_file, args.fusion_cluster_file)

print "test"
