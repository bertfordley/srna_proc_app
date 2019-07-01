#!/usr/bin/env python
"""
by: Robert C. Moseley
"""

import os
import json
import fuzzywuzzy.process
import subprocess
import argparse
import sys
import itertools
import pandas as pd


def trim_reads(fastq_in, fastq_trim_out, kit):

    # python 3 - can ran in parallel
    # cmd = "cutadapt -m 15 -u 3 -a {} -j 6 --quiet {} > {}".format(
    #     "A{10}", fastq_in, fastq_trim_out)
    # subprocess.call(cmd, shell=True)

    if kit.lower() == "clonetech":
        # python 2 - can't
        cmd = "cutadapt -m 15 -u 3 -a {} {} > {}".format("A{10}", fastq_in, fastq_trim_out)
        print(cmd)
        subprocess.call(cmd, shell=True)

    elif kit.lower() == "nextflex":

        outdir = os.path.split(fastq_trim_out)[0]
        trim_adapt_file = os.path.split(fastq_trim_out)[1][:-6] + "_adapt.fastq"
        trim_adapt_file = os.path.join(outdir, trim_adapt_file)
        # remove 3' adaptor "-a"  --  adapter not in reads but exact sequence is found in handle sequence
        cmd1 = "cutadapt -a TGGAATTCTCGGGTGCCAAGG -m 23 {} > {}".format(fastq_in, trim_adapt_file)
        print(cmd1)
        subprocess.call(cmd1, shell=True)

        # remove 4 bases from each end of read
        cmd2 = "cutadapt -u 4 -u -4 {} > {}".format(trim_adapt_file, fastq_trim_out)
        print()
        print(cmd2)
        subprocess.call(cmd2, shell=True)


def align_reads(fastq_file, index_file, sam_file, method, sensitivity):

    # user defined # of threads?
    cmd = "bowtie2 -p 6 --phred33 {} {} -x {} -U {} -S {}".format(
        method, sensitivity, index_file, fastq_file, sam_file)
    print()
    print(cmd)
    subprocess.call(cmd, shell=True)


def make_bam(sam_file):

    bam_file = sam_file[:-3] + "sorted.bam"
    cmd = "samtools view -b -F 0x100 -S {} | samtools sort -o {}".format(
        sam_file, bam_file)
    print()
    print(cmd)
    subprocess.call(cmd, shell=True)


def filter_spike_reads(bam_file, fastq_file, trimmed_file, seq_type, ids_file):

    if seq_type == "spike-ins":
        ids_file = fastq_file + "_" + seq_type + "_filtered_read_ids.txt"
        cmd1 = "samtools view -q 1 {} | cut -f1 > {}".format(
            bam_file, ids_file)
        print()
        print(cmd1)
        subprocess.call(cmd1, shell=True)

        filtered_fastq = trimmed_file[:-6] + "." + seq_type + "_filtered.fastq"
        cmd2 = "filterbyname.sh in={} names={} out={}".format(
            trimmed_file, ids_file, filtered_fastq)
        print()
        print(cmd2)
        subprocess.call(cmd2, shell=True)

    if seq_type == "grna-target-seqs":
        # bam_file not needed for gRNAs
        filtered_fastq = trimmed_file[:-6] + "." + seq_type + "_filtered.fastq"
        cmd2 = "filterbyname.sh in={} names={} out={}".format(
            trimmed_file, ids_file, filtered_fastq)
        print()
        print(cmd2)
        subprocess.call(cmd2, shell=True)


def merge_bam_files(fastq_file, bam_file1, bam_file2, bam_file3):

    merged_bam_files = fastq_file[:-9] + ".merged.sorted.bam"
    cmd = "samtools merge {} {} {} {}".format(
        merged_bam_files, bam_file1, bam_file2, bam_file3)
    print()
    print(cmd)
    subprocess.call(cmd, shell=True)


def match_parts(seq, gate_part_dict, grna_seq_list, thres):
    # match sgRNAs in the cicuit to reads in the fastq R1 and R2 files
    # fuzzy matching is used with user picking a threshold for percent
    # match.

    match = fuzzywuzzy.process.extractOne(seq, grna_seq_list)
    if match[1] >= thres:
        return (gate_part_dict[match[0]])


def read_fastq_seqs(filepath, line_type="seq"):
    # unzip and read fastq file

    if line_type == "header":
        with open(filepath, "rt") as fh:
            for seq_header, seq, qual_header, qual in itertools.zip_longest(*[fh] * 4):
                if len(seq.rstrip('\n')) >= 20:

                    yield seq_header.rstrip('\n')

    elif line_type == "seq":
        with open(filepath, "rt") as fh:
            for seq_header, seq, qual_header, qual in itertools.zip_longest(*[fh] * 4):
                if len(seq.rstrip('\n')) >= 20:

                    yield seq.rstrip('\n')


def count_grna_targets(filepath, gate_seqs_dict, threshold):
    grna_header_id_list = []
    grna_count_dic = {"r1": 0, "r2": 0, "r3": 0, "r5": 0, "r6": 0, "r7": 0, "r9": 0, "r10": 0}
    grna_seq_list = list(gate_seqs_dict.keys())

    header_gen = read_fastq_seqs(filepath, "header")
    seq_gen = read_fastq_seqs(filepath)
    print()
    print("fuzzywuzzy gRNA counter - file: {} threshold: {}".format(filepath, threshold))
    for (header, seq) in zip(header_gen, seq_gen):
        match = match_parts(seq, gate_seqs_dict, grna_seq_list, threshold)
        if match:
            grna_header_id_list.append(header.split(" ")[0][1:])
            grna_count_dic[match] += 1
    # print(grna_count_dic)
    return grna_header_id_list, grna_count_dic


if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog="small RNAseq fastq processing",
                                     description="Trims and aligns small RNAseq reads containing spike-ins")

    parser.add_argument("-fq", "--fastq_file",
                        type=str,
                        default=None,
                        required=True,
                        help="Must be forward read fastq file, including path")

    parser.add_argument("-si", "--spike_index_dir",
                        type=str,
                        default="indices/spike_ins/spike_ins",
                        required=True,
                        help="dir containing spike-in index files to align reads to")

    parser.add_argument("-ri", "--rna_index_dir",
                        type=str,
                        default="indices/srna_seqs/srna_seqs",
                        required=True,
                        help="dir containing srna sequences index files to align reads to")

    parser.add_argument("-t", "--threshold",
                        type=int,
                        default=90,
                        required=True,
                        help="threshold for fuzzy string matching reads to gRNA target sequences")

    parser.add_argument("-k", "--kit_tech",
                        type=str,
                        default="clonetech",
                        required=True,
                        help="small RNAseq kit used")

    parser.add_argument("-o", "--output_dir",
                        type=str,
                        default=None,
                        required=True,
                        help="output directory to write results")

    args = parser.parse_args()

    fastq_file = args.fastq_file
    fq_base_name = os.path.basename(fastq_file)
    spike_index_dir = args.spike_index_dir
    rna_index_dir = args.rna_index_dir
    threshold = args.threshold
    kit_tech = args.kit_tech
    output_dir = args.output_dir

    # arg check
    if not fastq_file:
        print("No read file given")
        sys.exit()
    if not spike_index_dir:
        print("No spike-in index directory given")
        sys.exit()
    if not rna_index_dir:
        print("No rna index directory given")
        sys.exit()
    if not kit_tech:
        print("Library prep kit not given")
        sys.exit()

    # make output directory
    os.mkdir(output_dir)


    # SPIKE-IN SEQUENCES#
    # ------------------------------------------------------------ #
    if kit_tech.lower() == "clonetech":
        # output path?
        trim_file = os.path.join(output_dir, fq_base_name[:-9] + ".trim.fastq")

        # use cutadapt to trim reads
        trim_reads(fastq_file, trim_file, kit_tech)
        spike_sam = os.path.join(output_dir, fq_base_name[:-9] + ".trim.spike-ins.sam")

    elif kit_tech.lower() == "nextflex":
        # output path?
        trim_file = os.path.join(output_dir, fq_base_name[:-9] + ".trim.fastq")

        # use cutadapt to trim reads
        trim_reads(fastq_file, trim_file, kit_tech)
        spike_sam = os.path.join(output_dir, fq_base_name[:-9] + ".trim.spike-ins.sam")

    # align spike-in reads using bowtie2
    align_reads(trim_file, spike_index_dir, spike_sam, "--end-to-end", "--very-sensitive")

    # make bam file for spike-ins
    make_bam(spike_sam)
    spike_bam = spike_sam[:-3] + "sorted.bam"

    # remove spike-in reads from fastq file
    filter_spike_reads(spike_bam, os.path.join(output_dir, fq_base_name[:-9]), trim_file, "spike-ins", False)
    # fastq file with no spike-in reads
    spike_filtered_fastq = trim_file[:-6] + ".spike-ins_filtered.fastq"


    # gRNA TARGET SEQUENCES #
    # ------------------------------------------------------------ #
    gate_part_seqs_dic = {'GGAACGTGATTGAATAACTT': "r1",
                          'ACCAACGCAAAAAGATTTAG': "r2",
                          'CATTGCCATACACCTTGAGG': "r3",
                          'GAAGTCAGTTGACAGAGTCG': "r5",
                          'GTGGTAACTTGCTCCATGTC': "r6",
                          'CTTTACGTATAGGTTTAGAG': "r7",
                          'GCAACCCACAAATATCCAGT': "r9",
                          'GTGACATAAACATTCGACTC': "r10"}

    grna_header_list, grna_counts_dic = count_grna_targets(spike_filtered_fastq, gate_part_seqs_dic, threshold)

    grna_filtered_ids_file = os.path.join(output_dir, fq_base_name[:-9] + "_grna-target-seqs_filtered_read_ids.txt")

    with open(grna_filtered_ids_file, 'w') as f:
        for header in grna_header_list:
            f.write("{}\n".format(header))

    output_df = pd.read_csv("indices/srna_grna_gene_lengths_kb.csv", index_col=0)

    # remove gRNA reads from fastq file
    filter_spike_reads(False, os.path.join(output_dir, fq_base_name[:-9]), spike_filtered_fastq,
                       "grna-target-seqs", grna_filtered_ids_file)
    # fastq file with no spike-in or gRNA reads
    grna_spike_filtered_fastq = spike_filtered_fastq[:-6] + ".grna-target-seqs_filtered.fastq"

    grna_json = json.dumps(grna_counts_dic)
    grna_json_file = open(os.path.join(output_dir, fq_base_name[:-9] + "_grna-counts.json"), "w")
    grna_json_file.write(grna_json)
    grna_json_file.close()

    # small RNA SEQUENCES #
    # ------------------------------------------------------------ #
    # make same for sRNA sequences
    # rna_sam = os.path.join(output_dir, fq_base_name[:-9] + ".trim.rna.sam")

    # align sRNA reads using bowtie2
    # align_reads(grna_spike_filtered_fastq, rna_index_dir, rna_sam, "--end-to-end", "--very-sensitive")

    # make bam for sRNA reads
    # make_bam(rna_sam)
    # rna_bam = rna_sam[:-3] + "sorted.bam"

    

    # merge spike-in, gRNA, and sRNA bam files
    # merge_bam_files(os.path.join(output_dir, fq_base_name), spike_bam, grna_bam, rna_bam)

    # get srna counts
