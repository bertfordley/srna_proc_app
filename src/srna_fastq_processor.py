#!/usr/bin/env python
"""
by: Robert C. Moseley
"""

import os
import shutil
import sys
import datetime
import subprocess
import argparse
import sys


def trim_reads(fastq_in, fastq_trim_out):

    cmd = "cutadapt -m 15 -u 3 -a {} -j 6 --quiet {} > {}".format(
        "A{10}", fastq_in, fastq_trim_out)
    subprocess.call(cmd, shell=True)


def align_reads(fastq_file, index_file, sam_file):

    # user defined # of threads?
    cmd = "bowtie2 -p 6 --phred33 --very-sensitive -x {} -U {} -S {}".format(
        index_file, fastq_file, sam_file)
    subprocess.call(cmd, shell=True)


def make_bam(sam_file):

    bam_file = sam_file[:-3] + ".sorted.bam"
    cmd = "samtools view -b -F 0x100 -S {} | samtools sort -o {}".format(
        sam_file, bam_file)
    subprocess.call(cmd, shell=True)


def filter_spike_reads(bam_file, fastq_file, trimmed_file):

    ids_file = fastq_file[:"-9"] + "filtered_read_ids.txt"
    cmd1 = "samtools view -q 1 {} | cut -f1 > {}".format(
        bam_file, ids_file)
    subprocess.call(cmd1, shell=True)

    filtered_fastq = trimmed_file[:-6] + "_filtered.fastq"
    cmd2 = "filterbyname.sh in={} names={} out={}".format(
        trimmed_file, ids_file, filtered_fastq)
    subprocess.call(cmd2, shell=True)


def merge_bam_files(fastq_file, bam_file1, bam_file2):

    merged_bam_files = fastq_file[:-9] + ".merged.sorted.bam"
    cmd = "samtools merge {} {} {}".format(
        merged_bam_files, bam_file1, bam_file2)
    subprocess.call(cmd, shell=True)


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
                        default=None,
                        required=True,
                        help="dir containing spike-in index files to align reads to")

    parser.add_argument("-ri", "--rna_index_dir",
                        type=str,
                        default=None,
                        required=True,
                        help="dir containing srna and grna index files to align reads to")

    parser.add_argument("-o", "--output_dir",
                        type=str,
                        default=None,
                        required=True,
                        help="output directory to write results")

    args = parser.parse_args()

    fastq_file = args.fastq_file
    spike_index_dir = args.spike_index_dir
    rna_index_dir = args.rna_index_dir
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

    ## Spike-ins

    # output path?
    trim_file = fastq_file[:-9] + "trim.fastq"

    # use cutadapt to trim reads
    trim_reads(fastq_file, trim_file)
    spike_sam = fastq_file[:-9] + "trim.spike-ins.sam"

    # align spike-in reads using bowtie2
    align_reads(trim_file, spike_index_dir, spike_sam)

    # make bam file for spike-ins
    make_bam(spike_sam)
    spike_bam = spike_sam[:-3] + "sorted.bam"

    # remove spike-in reads from fastq file
    filter_spike_reads(spike_bam, fastq_file, trim_file)
    # fastq file with no spike-in reads
    filtered_fastq = trim_file[:-6] + "_filtered.fastq"

    ## gRNAs and sRNAs

    rna_sam = fastq_file[:-9] + "_rna.sam"
    # align gRNA/sRNA reads using bowtie2
    align_reads(filtered_fastq, rna_index_dir, rna_sam)

    # make bam for gRNA/sRNA
    make_bam(rna_sam)
    rna_bam = rna_sam[:-3] + "sorted.bam"

    # merge spike-in and gRNA/sRNA bam files
    merge_bam_files(fastq_file, spike_bam, rna_bam)
