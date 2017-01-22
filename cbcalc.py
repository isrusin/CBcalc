#! /usr/bin/python

import argparse as ap
import sys
from os.path import basename

import cbclib.sites
from cbclib.counts import Counts

__version__ = "1.1"

def make_output_stubs(methods):
    """Make headers and row stub for output table.

    Arguments:
        methods -- a list of method abbreviations from
                  sbslib.sites.wrappers

    Returns:
        headers -- output table headers as string
        row_stub -- output table row stub to use with str.format()
    """
    headers = "Sequence ID\tSite\tObserved\t"
    row_stub = "{id}\t{site}\t{num:d}\t"
    for method in methods:
        headers += "Expected (%s)\tRatio (%s)\t" % (method, method)
        row_stub += "{%sexp:.2f}\t{%srat:.3f}\t" % (method, method)
    headers += "Total\n"
    row_stub += "{total:.0f}\n"
    return headers, row_stub

def wrap_sites(raw_sites, methods, maxlen=10):
    """Wrap raw sites with wrappers implementing specified methods.

    Arguments:
        raw_sites -- a list of sites as strings
        methods -- a list of method abbreviations from
                   cbclib.sites.wrappers
        maxlen (optional) -- site length cutoff, default 10

    Returns:
        wrapped_sites -- a list of tuples, each contains wrapped
                         versions of a site for each method
        unwrapped_sites -- a dict which contains reason to skip for every
                           unwrapped raw site
    """
    wrapped_sites = []
    unwrapped_sites = {}
    for raw_site in raw_sites:
        wrapped_site = []
        for method in methods:
            Wrapper = cbclib.sites.get_wrapper_by_abbr(method)
            try:
                wrapped_site.append(Wrapper(raw_site, maxlen))
            except ValueError as error:
                unwrapped_sites[raw_site] = error.message
                break
        else:
            wrapped_sites.append(tuple(wrapped_site))
    return wrapped_sites, unwrapped_sites

def cbcalc(sid, row_stub, sites, counts, methods):
    """Calculate values for sites and return formatted output table rows.

    Arguments:
        sid -- sequence ID to put into 'ID' column of the output table
        row_stub -- stub of table row to use in str.format
        sites -- a list of wrapped sites obtained with wrap_sites()
        counts -- word counts calculated with cbclib.count.calc_all()
        methods -- a list of method abbreviations from
                   cbclib.sites.wrappers

    Returns:
        rows -- a list of formatted rows of the output table
    """
    rows = []
    vals = {"id": sid}
    for wrapped in sites:
        index = 0
        wsite = wrapped[index]
        vals["site"] = wsite.str_init
        vals["total"] = counts.get_total(wsite.struct_hash)
        obs = wsite.calc_observed(counts)
        vals["num"] = obs
        for method in methods:
            wsite = wrapped[index]
            exp = wsite.calc_expected(counts)
            vals["%sexp" % method] = exp
            vals["%srat" % method] = obs / (exp or float("NaN"))
            index += 1
        rows.append(row_stub.format(**vals))
    return rows

def main(argv=None):
    parser = ap.ArgumentParser(
        description="CBcalc - Compositional Bias calculation.",
        add_help=False, usage="\n    %(prog)s ".join([
            "", "--help", "--version", "[-soMPK] FASTA [FASTA ...]",
            "[-soMPK] PATH -i LIST"
        ])
    )
    parser.add_argument(
        "inseq", metavar="FASTA", nargs="+",
        help="Input .fasta file, may be gzipped."
    )
    parser.add_argument(
        "-i", "--id", dest="insids", metavar="LIST", type=ap.FileType("r"),
        help="""Input LIST file with sequence IDs; makes CBcalc to treat
        the first (and the only in the case) positional argument as PATH
        stub."""
    )
    parser.add_argument(
        "_none", metavar="PATH", nargs="?",# only to make up help message
        help="""Input .fasta files path stub with '{}' placeholder for
        sequence IDs, which are listed in a LIST file specified with
        -i/--id option."""
    )
    parser.add_argument(
        "-s", "--site", dest="instl", metavar="LIST",
        type=ap.FileType('r'), default=sys.stdin,
        help="File with a list of sites, one-per-line, default is stdin."
    )
    parser.add_argument(
        "-o", "--out", dest="outsv", metavar="FILE",
        type=ap.FileType('w'), default=sys.stdout,
        help="Output tabular (.tsv) file, default is stdout."
    )
    method_group = parser.add_argument_group(
        description="""Folowing arguments allow
        to select methods of expected frequency calculation. Order of
        the arguments determines column order of the output file. If
        no method is specified, default set (mmax, pevzner, karlin)
        will be used."""
    )
    method_group.add_argument(
        "-B", "--bernoulli", dest="methods", action="append_const",
        const="B", help="Bernoulli model based method"
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
    parser.add_argument(
        "-v", "--version", action="version",
        version="%(prog)s " + __version__, help=ap.SUPPRESS
    )
    parser.add_argument(
        "-h", "-?", "-help", "--help", action="help",
        help=ap.SUPPRESS
    )
    args = parser.parse_args(argv)
    if args.insids is not None:
        with args.insids as insids:
            sids = sorted(set(insids.read().split()))
        seqs = [args.inseq[0].format(sid) for sid in sids]
    else:
        seqs = args.inseq
        sids = [basename(seq).split(".fasta")[0] for seq in seqs]
    methods = args.methods or ["M", "P", "K"]
    headers, row_stub = make_output_stubs(methods)
    methods = sorted(set(methods))
    with args.instl as instl:
        raw_sites = instl.read().split()
    sites, unwrapped = wrap_sites(raw_sites, methods)
    structs = cbclib.sites.get_structs([_s[0] for _s in sites], methods)
    with args.outsv as outsv:
        outsv.write(headers)
        for sid, seq in zip(sids, seqs):
            counts = Counts(seq, structs)
            rows = cbcalc(sid, row_stub, sites, counts, methods)
            outsv.writelines(rows)


if __name__ == "__main__":
    sys.exit(main())
