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

def cbcalc(outsv, ouline, sites, counts, methods):
    for wrapped in sites:
        vals = dict()
        wsite = wrapped[0]
        vals["Site"] = wsite.str_init
        vals["Total"] = counts.get_total(wsite.struct)
        obs = wsite.calc_observed(counts)
        vals["No"] = obs
        for method in methods:
            wsite = wrapped.pop(0)
            exp = wsite.calc_expected(counts)
            vals["%ce" % method] = exp
            vals["%cr" % method] = obs / (exp or float("NaN"))
        outsv.write(ouline.format(**vals))

if __name__ == "__main__":
    parser = ap.ArgumentParser(
            description="Contrast calculation", usage="cbcalc.py [-h] " +
            "FASTA [-s file] [-o file] [method(s)]"
            )
    parser.add_argument(
            "inseq", metavar="FASTA",
            help="Input .fasta file, may be gzipped."
            )
    parser.add_argument(
            "-s", "--sites", dest="instl", metavar="file",
            type=ap.FileType('r'), default=sys.stdin,
            help="Input list of sites, one-per-line."
            )
    parser.add_argument(
            "-o", "--out", dest="outsv", metavar="file",
            type=ap.FileType('w'), default=sys.stdout,
            help="Output tabular (.tsv) file, default is stdout"
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
    args = parser.parse_args()
    methods = args.methods or ["M", "P", "K"]
    title = "Site\tNo\t"
    ouline = "{Site}\t{No:d}\t"
    for method in methods:
        title += "%se\t%sr\t" % (method, method)
        ouline += "{%se:.2f}\t{%sr:.3f}\t" % (method, method)
    title += "Total\n"
    ouline += "{Total:.0f}\n"
    methods = sorted(set(methods))
    sites, structs = load_sites(args.instl, methods)
    with cnt.calc_all(args.inseq, structs) as counts, args.outsv as outsv:
        outsv.write(title)
        cbcalc(outsv, ouline, sites, counts, methods)

