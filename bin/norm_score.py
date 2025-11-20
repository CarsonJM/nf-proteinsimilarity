#!/usr/bin/env python

import argparse
import csv
import gzip
import multiprocessing as mp
import os
import psutil
import shutil
import signal
import subprocess as sp
import time

from collections import defaultdict

def parse_args(args=None):
    description = "Calculate the normalized bitscore for a genome."
    epilog = "Example usage: python uhvdb_normscore.py --help"

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-i",
        "--input",
        help="Path to TSV created by DIAMOND.",
    )
    parser.add_argument(
        "-s",
        "--self_score",
        help="Path to self score TSV created using self_score.py.",
    )
    parser.add_argument(
        "-t",
        "--threads",
        help="Number of threads to use for parallel processing.",
        type=int
    )
    parser.add_argument(
        "-m",
        "--min_score",
        help="Minimum normscore to output.",
        type=int
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output TSV containing normalized bitscore.",
    )
    parser.add_argument('--version', action='version', version='1.0.0')
    return parser.parse_args(args)


def mean(values):
    return sum(values) / len(values)


def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def terminate_tree(pid, including_parent=True):
    parent = psutil.Process(pid)
    for child in parent.children(recursive=True):
        child.terminate()
    if including_parent:
        parent.terminate()


def run_shell(cmd):
    p = sp.Popen(cmd, shell=True)
    return p.wait()

        
def parallel(function, arguments_list, threads):
    pool = mp.Pool(threads, init_worker)
    try:
        results = []
        for arguments in arguments_list:
            p = pool.apply_async(function, args=arguments)
            results.append(p)
        pool.close()
        while True:
            if all(r.ready() for r in results):
                return [r.get() for r in results]
            time.sleep(1)
    except KeyboardInterrupt:
        pid = os.getpid()
        terminate_tree(pid)


def split_dmnd(inpath, outdir, num_splits, ext=""):

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    total_size = os.stat(inpath).st_size
    split_size = int(total_size / num_splits)

    last_id = None
    split_num = 1
    cursize = 0
    outfile = open(os.path.join(outdir, str(split_num)), "w")
    if inpath.endswith('.gz'):
        infile = gzip.open(inpath, 'rt')
    else:
        infile = open(inpath, 'r')
    # for line in infile:
    #     outfile.write(line)
    #     cursize += len(line)
    #     last_id = line.split()[0].rsplit("_", 1)[0]
    for line in infile:
        cur_id = line.split()[0].rsplit("_", 1)[0]
        if cursize > split_size and cur_id != last_id:
            split_num += 1
            print(f"Created split {split_num}", flush=True)
            cursize = 0
            outfile = open(os.path.join(outdir, str(split_num)), "w")
        outfile.write(line)
        cursize += len(line)
        last_id = cur_id
    outfile.close()
    infile.close()


def yield_diamond_hits(diamond):
    with open(diamond) as f:
        try:
            hits = [next(f).split()]
        except StopIteration:
            return
        for line in f:
            r = line.split()
            query = r[0].rsplit("_", 1)[0]
            last = hits[-1][0].rsplit("_", 1)[0]
            if query != last:
                yield last, hits
                hits = []
            hits.append(r)
        if len(hits) > 0:
            last = hits[-1][0].rsplit("_", 1)[0]
            yield last, hits


def split_hits(hits):
    target_to_hits = defaultdict(list)
    for hit in hits:
        tname = hit[1].rsplit("_", 1)[0]
        target_to_hits[tname].append(hit)
    return target_to_hits


def best_blast_hits(hits, query_key=0, score_key=-1):
    bhits = {}
    for hit in hits:
        if hit[query_key] not in bhits:
            bhits[hit[query_key]] = hit
        elif float(hit[score_key]) > float(bhits[hit[query_key]][score_key]):
            bhits[hit[query_key]] = hit
    return list(bhits.values())


def aai_main(inpath, outpath, selfpath, min_score):
    selfaai = {}
    if selfpath.endswith('.gz'):
        selfpath_read = gzip.open(selfpath, 'rt')
    else:
        selfpath_read = open(selfpath, 'r')
    for r in csv.DictReader(selfpath_read, delimiter="\t"):
        selfaai[r["genome_id"]] = float(r["selfscore"])
    with open(outpath, "w") as out:
        header = ["query", "reference", "hits", "aai", "raw_score", "norm_score"]
        out.write("\t".join(header) + "\n")
        for qname, hits in yield_diamond_hits(inpath):
            target_to_hits = split_hits(hits)
            for tname, thits in target_to_hits.items():
                if qname == tname:
                    continue
                bhits = best_blast_hits(thits)
                aai = mean([float(_[2]) for _ in bhits])
                score = sum([float(_[-1]) for _ in bhits])
                norm = 100 * score / selfaai[qname]
                if norm < min_score:
                    continue
                row = [
                    qname,
                    tname,
                    len(bhits),
                    round(aai, 2),
                    round(score, 2),
                    round(norm, 2),
                ]
                out.write("\t".join([str(_) for _ in row]) + "\n")

def main(args=None):
    args = parse_args(args)

    split_dmnd(args.input, f"./{args.input}_dmnddir/", args.threads)

    if not os.path.exists(f"./{args.input}_aaidir/"):
        os.makedirs(f"./{args.input}_aaidir/")
            
    argument_list = []
    for file in os.listdir(f"./{args.input}_dmnddir/"):
        inpath = os.path.join(f"./{args.input}_dmnddir/", file)
        outpath = os.path.join(f"./{args.input}_aaidir/", file)
        argument_list.append([inpath, outpath, args.self_score, args.min_score])
    parallel(aai_main, argument_list, threads=args.threads)

    with open(args.output, "w") as out:

        for file in os.listdir(f"./{args.input}_aaidir/"):
            handle = open(os.path.join(f"./{args.input}_aaidir/", file))
            out.write(next(handle))
            break

        for file in os.listdir(f"./{args.input}_aaidir/"):
            handle = open(os.path.join(f"./{args.input}_aaidir/", file))
            next(handle)
            for line in handle:
                out.write(line)
            handle.close()

    shutil.rmtree(f"./{args.input}_dmnddir/")
    shutil.rmtree(f"./{args.input}_aaidir/")

if __name__ == "__main__":
    main()
