#! /usr/bin/python

import argparse as ap
import sys

import cbclib.counts as cnt
import cbclib.sites as st


def load_sites(instl, flags, len_cutoff=10):
    sites = []
    for_structs = []
    wrappers = zip(flags, [st.MarkovSite, st.PevznerSite, st.KarlinSite])
    with instl:
        for line in instl:
            site = line.strip("\n\r\t Nn-.")
            wrapped = []
            for flag, wrapper in wrappers:
                if flag:
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

if __name__ == "__main__":
    parser = ap.ArgumentParser(
            description="Contrast calculation", usage="usage: cbcalc.py " +
            "[-h] FASTA [-s file] [-o file] [method(s)]"
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
            const="mmax", help="Mmax based method"
            )
    method_group.add_argument(
            "-P", "--pevzner", dest="methods", action="append_const",
            const="pevzner", help="Pevzner's method"
            )
    method_group.add_argument(
            "-K", "--karlin", dest="methods", action="append_const",
            const="karlin", help="Karlin's method"
            )
    args = parser.parse_args()
    methods_list = ["mmax", "pevzner", "karlin"]
    methods = args.methods
    if not methods:
        methods = methods_list
    flags = []
    for method in methods_list:
        flags.append(method in methods)
    sites, structs = load_sites(args.instl, flags)
    counts = cnt.calc_all(args.inseq, structs)
    title = "Site\tNo\t"
    ouline = "{Site}\t{No:d}\t"
    for method in methods:
        substitution = tuple(method[0].upper() * 2)
        title += "%ce\t%cr\t"  % substitution
        ouline += "{%ce:.2f}\t{%cr:.3f}\t" % substitution
    title += "Total number\n"
    ouline += "{Total:.0f}\n"
    with args.outsv as outsv:
        outsv.write(title)
        for wrapped in sites:
            vals = dict()
            wsite = wrapped[0]
            vals["Site"] = wsite.str_init
            vals["Total"] = counts.get_total(wsite.struct)
            vals["No"] = wsite.calc_observed(counts)
            if flags[0]: #Mmax
                wsite = wrapped.pop(0)
                vals["Me"] = wsite.calc_expected(counts)
                vals["Mr"] = vals["No"] / (vals["Me"] or float("NaN"))
            if flags[1]: #Pevzner
                wsite = wrapped.pop(0)
                vals["Pe"] = wsite.calc_expected(counts)
                vals["Pr"] = vals["No"] / (vals["Pe"] or float("NaN"))
            if flags[2]: #Karlin
                wsite = wrapped.pop(0)
                vals["Ke"] = wsite.calc_expected(counts)
                vals["Kr"] = vals["No"] / (vals["Ke"] or float("NaN"))
            outsv.write(ouline.format(**vals))

