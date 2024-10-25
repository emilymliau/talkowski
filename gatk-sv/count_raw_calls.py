"""
GATK-SV QC: Count Raw SV Calls
---
Author: Emily M. Liau
Last Updated: October 25, 2024

Description: This script aggregates the number of raw structural variant (SV) calls from the primary SV 
callers (Manta, Wham, MELT) to supplement the existing metrics in preliminary quality control (QC) phase 
that follows the 02-EvidenceQC stage of the GATK-SV pipeline.

Input(s):
- gzipped VCF files for each SV caller, generated as output from 01-GatherSampleEvidence
- files should be organized in a single parent directory, with separate subdirectories for each SV caller

Output(s):
- individual TSV files for each SV caller, saved within their respective subdirectories
"""

#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

# import packages
import os
import argparse
import gzip
from collections import defaultdict

# extracts chromosome information from VCF files
def chr_extract(pin, chr):
    out = ''
    for i in pin[7].split(';'):
        if i.split('=')[0] == chr:
            out = i.split('=')[1]
    return out

# reads a gzipped VCF file and counts the occurrences of different SV types
def vcf_stat_readin(filename):
    fin = gzip.open(filename, 'rb')
    sv_counts = defaultdict(int)  # dictionary to store counts of SV types
    sv_types = set()  # initialize set to identify unique SV types observed in VCF file
    for line in fin:
        pin = line.strip().decode().split()
        if not pin[0][0] == '#':
            svtype = chr_extract(pin, 'SVTYPE')
            if svtype:
                sv_counts[svtype] += 1  # track observed frequency of SV types
                sv_types.add(svtype)  # add all observed SV types to set
    fin.close()
    return sv_counts, sv_types

# writes SV counts for a sample to an output file
def write_output(output_name, sv_counts, sv_types, sample_name):
    # create file with SV types as column names
    if not os.path.isfile(output_name):
        with open(output_name, 'w') as fo:
            header = ['SAMPLE'] + sorted(sv_types)
            print('\t'.join(header), file=fo)
    
    # add new row for each sample
    with open(output_name, 'a') as fo:
        row = [sample_name]
        for sv_type in sorted(sv_types):
            row.append(str(sv_counts.get(sv_type, 0)))  # add SV counts for each SV type
        print('\t'.join(row), file=fo)

# sorts the TSV file by the 'SAMPLE' column and rewrites the information to the file
def sort_tsv(output_file):
    with open(output_file, 'r') as fo:
        lines = fo.readlines()

    # extract the header and the rest of the lines
    header = lines[0]
    data = lines[1:]

    # sort the data based on the first column (SAMPLE)
    sorted_data = sorted(data, key=lambda x: x.split('\t')[0])

    # write the sorted data back to the file
    with open(output_file, 'w') as fo:
        fo.write(header)
        fo.writelines(sorted_data)

# process all gzipped VCF files in specified directory
def process_directory(vcf_dir, output_file):
    # extract all .vcf.gz files in given directory
    vcf_files = [os.path.join(vcf_dir, f) for f in os.listdir(vcf_dir) if f.endswith('.vcf.gz')]

    # iterate through each VCF file in the specified directory
    for vcf_file in vcf_files:
        sample_name = os.path.basename(vcf_file).split('.')[0]
        sv_counts, sv_types = vcf_stat_readin(vcf_file)
        
        # add SV counts for the current sample to the output TSV file
        write_output(output_file, sv_counts, sv_types, sample_name)

    # sort the TSV file after processing all VCF files
    sort_tsv(output_file)

def main():
    parser = argparse.ArgumentParser("count_raw_calls.py")
    parser.add_argument("base_dir", type=str, help="parent directory containing SV caller-specific "
                        "subdirectories (manta_vcfs, melt_vcfs, wham_vcfs)")
    args = parser.parse_args()

    # define the subdirectories containing gzipped VCF files from each SV caller
    subdirectories = {
        "manta": "manta_vcfs",
        "melt": "melt_vcfs",
        "wham": "wham_vcfs"
    }

    # process each subdirectory of VCF files
    for algorithm, subdir in subdirectories.items():
        vcf_dir = os.path.join(args.base_dir, subdir)
        output_file = os.path.join(vcf_dir, f"{algorithm}_sv_counts.tsv")
        
        # process the current directory and write output to a separate TSV file for each SV caller
        process_directory(vcf_dir, output_file)


if __name__ == '__main__':
    main()
