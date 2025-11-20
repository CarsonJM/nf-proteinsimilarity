#!/usr/bin/env python

import argparse
import csv
import gzip

def parse_args(args=None):
    description = "Calculate AAI self-alignment score for a genome."
    epilog = "Example usage: python self_score.py --help"

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-i",
        "--input",
        help="Path to TSV created by DIAMOND (should include self-alignments).",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output TSV containing AAI self-alignment score.",
    )
    parser.add_argument('--version', action='version', version='1.0.0')
    return parser.parse_args(args)


def calculate_aai(diamond_file):
    data = {}

    # initialize dict
    # only use gzip if needed
    if diamond_file.endswith('.gz'):
        diamond_read = gzip.open(diamond_file, 'rt')
    else:
        diamond_read = open(diamond_file, 'r')

    line_num = 0

    for r in csv.reader(diamond_read, delimiter="\t"):
        genome = r[0].rsplit("_", 1)[0]
        data[genome] = {"genes": 0, "selfscore": 0}
        line_num += 1
        if line_num % 1000000 == 0:
            print(f"Processed {line_num} lines...", flush=True)

    # reset file handle
    diamond_read.close()


    # only use gzip if needed
    if diamond_file.endswith('.gz'):
        diamond_read = gzip.open(diamond_file, 'rt')
    else:
        diamond_read = open(diamond_file, 'r')

    line_num = 0

    for r in csv.reader(diamond_read, delimiter="\t"):
        genome = r[0].rsplit("_", 1)[0]
        line_num += 1
        # print(r)
        if r[0] == r[1]:
            data[genome]["selfscore"] += float(r[-1])
            data[genome]["genes"] += 1

        if line_num % 1000000 == 0:
            print(f"Processed {line_num} lines...", flush=True)

    return data

def main(args=None):
    args = parse_args(args)

    # Modify the file
    data = calculate_aai(args.input)

    with open(args.output, "w") as out:
        header = ["genome_id", "genes", "selfscore"]
        out.write("\t".join(header) + "\n")
        for id in data:
            # print(data)
            rec = [id, data[id]["genes"], data[id]["selfscore"]]
            out.write("\t".join([str(_) for _ in rec]) + "\n")

if __name__ == "__main__":
    main()
