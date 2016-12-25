#! /usr/bin/python

import argparse as ap
import sys
from os import basename

import cbclib.counts
import cbclib.sites


def make_output_stubs(methods):
    """Make headers and row stub for output table.

    Arguments:
        methods -- a list of method abbreviations from
                  sbslib.sites.wrappers

    Returns:
        headers -- output table headers as string
        row_stub -- output table row stub to use with str.format()
    """
    headers = "ID\tSite\tNo\t"
    row_stub = "{id}\t{site}\t{num:d}\t"
    for method in methods:
        headers += "%se\t%sr\t" % (method, method)
        row_stub += "{%se:.2f}\t{%sr:.3f}\t" % (method, method)
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
            Wrapper = cbclib.sites.wrappers[method]
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
        vals["total"] = counts.get_total(wsite.struct)
        obs = wsite.calc_observed(counts)
        vals["num"] = obs
        for method in methods:
            wsite = wrapped[index]
            exp = wsite.calc_expected(counts)
            vals["%ce" % method] = exp
            vals["%cr" % method] = obs / (exp or float("NaN"))
            index += 1
        rows.append(row_stub.format(**vals))
    return rows

def main(argv=None):
    parser = ap.ArgumentParser(description="Contrast calculation")
    parser.add_argument(
            "-s", "--sites", dest="instl", metavar="file",
            type=ap.FileType('r'), default=sys.stdin,
            help="Input list of sites, one-per-line, default is stdin."
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
            file as sequence ID (basenames without .fasta* extension)."""
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
            help="""Input file with a list of sequence IDs, IDs should not
            contain any whitespace symbols."""
            )
    args = parser.parse_args(argv)
    ispath = "sids" in vars(args)
    if ispath:
        with open(args.sids) as insids:
            sids = sorted(set(insids.read().split()))
        seqs = [args.inseq.format(sid) for sid in sids]
    else:
        seqs = args.inseq
        sids = [basename(seq).split(".fasta")[0] for seq in seqs]
    methods = args.methods or ["M", "P", "K"]
    headers, row_stub = make_output_stubs(methods)
    methods = sorted(set(methods))
    with args.instl as instl:
        raw_sites = instl.read().split()
    sites, unwrapped = wrap_sites(raw_sites, methods)
    structs = cbclib.sites.get_structs(s for w in sites for s in w)
    with args.outsv as outsv:
        outsv.write(headers)
        for sid, seq in zip(sids, seqs):
            with cbclib.counts.calc_all(seq, structs) as counts:
                rows = cbcalc(sid, row_stub, sites, counts, methods)
            outsv.writelines(rows)


if __name__ == "__main__":
    sys.exit(main())
