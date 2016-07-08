#! /usr/bin/python

import argparse as ap
import sys

import cbclib.counts as cnt
import cbclib.sites as st


def load_sites(instl, methods, len_cutoff=10):
    sites = []
    for_structs = []
    with instl:
        for line in instl:
            site = line.strip("\n\r\t Nn-.")
            wrapped = []
            for method in methods:
                wrapper = st.wrappers[method]
                try:
                    wsite = wrapper(site)
                    length = wsite.L
                    if length > len_cutoff:
                        break
                    wrapped.append(wsite)
                except ValueError:
                    break
            else:
                sites.append(wrapped)
                for_structs.extend(wrapped)
                continue
            sys.stderr.write("%s is too long, skipped.\n" % site)
    structs = st.get_structs(for_structs)
    return sites, structs

def cbcalc(sid, outsv, ouline, sites, counts, methods):
    vals = {"ID": sid}
    index = 0
    for wrapped in sites:
        wsite = wrapped[index]
        vals["Site"] = wsite.str_init
        vals["Total"] = counts.get_total(wsite.struct)
        obs = wsite.calc_observed(counts)
        vals["No"] = obs
        for method in methods:
            wsite = wrapped[index]
            exp = wsite.calc_expected(counts)
            vals["%ce" % method] = exp
            vals["%cr" % method] = obs / (exp or float("NaN"))
            index += 1
        outsv.write(ouline.format(**vals))

if __name__ == "__main__":
    parser = ap.ArgumentParser(description="Contrast calculation")
    parser.add_argument(
            "-s", "--sites", dest="instl", metavar="file",
            type=ap.FileType('r'), default=sys.stdin,
            help="Input list of sites, one-per-line."
            )
    parser.add_argument(
            "-o", "--out", dest="outsv", metavar="file",
            type=ap.FileType('w'), default=sys.stdout,
            help="Output tabular (.tsv) file, default is stdout."
            )
    method_group = parser.add_argument_group(
            title="Method arguments", description="""Arguments that allow
            to select methods of expected frequency calculation. Order of
            the arguments determines column order of the output file. If
            no method is specified, default set (mmax, pevzner, karlin)
            will be used."""
            )
    method_group.add_argument(
            "-M", "--mmax", dest="methods", action="append_const",
            const="M", help="Mmax based method"
            )
    method_group.add_argument(
            "-P", "--pevzner", dest="methods", action="append_const",
            const="P", help="Pevzner's method"
            )
    method_group.add_argument(
            "-K", "--karlin", dest="methods", action="append_const",
            const="K", help="Karlin's method"
            )
    subparsers = parser.add_subparsers(
            title="Input mode", metavar="MODE",
            help="Either 'file' or 'path'."
            )
    file_parser = subparsers.add_parser(
            "file", help="""In 'file' mode you should specify input
            .fasta(.gz) files by name which will be used in the output
            file as sequence ID (without .fasta* extension)."""
            )
    file_parser.add_argument(
            "inseq", metavar="FASTA", nargs="+",
            help="Input .fasta file, may be gzipped."
            )
    path_parser = subparsers.add_parser(
            "path", help="""In 'path' mode you should specify a single
            path with {} placeholder for sequence ID which will be
            extracted from a file (-i/--id option) and will be used in the
            output file."""
            )
    path_parser.add_argument(
            "inseq", metavar="PATH", help="""Input .fasta files path,
            use {} as placeholder for sequence IDs, specified with
            -i/--id option."""
            )
    path_parser.add_argument(
            "-i", "--id", dest="sids", metavar="file", required=True,
            help="Input file with a list of sequence IDs, one-per-line."
            )
    args = parser.parse_args()
    methods = args.methods or ["M", "P", "K"]
    title = "ID\tSite\tNo\t"
    ouline = "{ID}\t{Site}\t{No:d}\t"
    for method in methods:
        title += "%se\t%sr\t" % (method, method)
        ouline += "{%se:.2f}\t{%sr:.3f}\t" % (method, method)
    title += "Total\n"
    ouline += "{Total:.0f}\n"
    methods = sorted(set(methods))
    sites, structs = load_sites(args.instl, methods)
    ispath = "sids" in vars(args)
    if ispath:
        with open(args.sids) as insids:
            sids = sorted(set(insids.read().strip().split("\n")))
        seqs = [args.inseq.format(sid) for sid in sids]
    else:
        seqs = args.inseq
        sids = [seq.split(".fasta")[0] for seq in seqs]
    with args.outsv as outsv:
        outsv.write(title)
        for sid, seq in zip(sids, seqs):
            with cnt.calc_all(seq, structs) as counts:
                cbcalc(sid, outsv, ouline, sites, counts, methods)

