#! /usr/bin/python

"""Compositional Bias calculation tool."""

import argparse as _ap
import sys as _sys
from os.path import basename as _basename

import cbclib.sites as _sites
from cbclib.counts import Counts as _Counts
from cbclib import __version__

def make_output_stubs(methods):
    """Make headers and row stub for output table.

    Arguments:
        methods -- a list of method abbreviations from cbclib.sites

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
        methods -- a list of method abbreviations from cbclib.sites
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
            Wrapper = _sites.get_wrapper_by_abbr(method)
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
        counts -- cbclib.counts.Counts instance
        methods -- a list of method abbreviations from cbclib.sites

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
    """Main function.

    Parse CLI arguments and execute `cbcalc` function.
    """
    parser = _ap.ArgumentParser(
        description="CBcalc - Compositional Bias calculation.",
        add_help=False, usage="\n    %(prog)s ".join([
            "", "--help", "--version", "[-soBMPK] FASTA [FASTA ...]",
            "[-soBMPK] PATH -i LIST"
        ])
    )
    parser.add_argument(
        "inseq", metavar="FASTA", nargs="+",
        help="Input .fasta file, may be gzipped."
    )
    parser.add_argument(
        "-i", "--id", dest="insids", metavar="LIST",
        type=_ap.FileType("r"), help="""Input LIST file with sequence IDs;
        makes CBcalc to treat the first (and the only in the case)
        positional argument as PATH stub."""
    )
    parser.add_argument(
        "_none", metavar="PATH", nargs="?",# only to make up help message
        help="""Input .fasta files path stub with '{}' placeholder for
        sequence IDs, which are listed in a LIST file specified with
        -i/--id option."""
    )
    parser.add_argument(
        "-s", "--site", dest="instl", metavar="LIST",
        type=_ap.FileType('r'), default=_sys.stdin,
        help="File with a list of sites, one-per-line, default is stdin."
    )
    parser.add_argument(
        "-o", "--out", dest="outsv", metavar="FILE",
        type=_ap.FileType('w'), default=_sys.stdout,
        help="Output tabular (.tsv) file, default is stdout."
    )
    method_group = parser.add_argument_group(
        description="""Folowing arguments determine which methods of
        expected frequency calculation to use. The order of the
        arguments defines the column order of the output file. If
        no method is specified, default set (mmax, pevzner, karlin)
        will be used."""
    )
    method_group.add_argument(
        "-B", "--bernoulli", dest="methods", action="append_const",
        const=_sites.Site.abbr, help="Bernoulli model based method"
    )
    method_group.add_argument(
        "-M", "--mmax", dest="methods", action="append_const",
        const=_sites.MarkovSite.abbr, help="Mmax based method"
    )
    method_group.add_argument(
        "-P", "--pevzner", dest="methods", action="append_const",
        const=_sites.PevznerSite.abbr, help="Pevzner's method"
    )
    method_group.add_argument(
        "-K", "--karlin", dest="methods", action="append_const",
        const=_sites.KarlinSite.abbr, help="Karlin's method"
    )
    parser.add_argument(
        "-v", "--version", action="version",
        version=("CBcalc v%s" % __version__), help=_ap.SUPPRESS
    )
    parser.add_argument(
        "-h", "-?", "-help", "--help", action="help",
        help=_ap.SUPPRESS
    )
    args = parser.parse_args(argv)
    if args.insids is not None:
        with args.insids as insids:
            sids = sorted(set(insids.read().split()))
        seqs = [args.inseq[0].format(sid) for sid in sids]
    else:
        seqs = args.inseq
        sids = [_basename(seq).split(".fasta")[0] for seq in seqs]
    methods = args.methods or [
        _sites.MarkovSite.abbr,
        _sites.PevznerSite.abbr,
        _sites.KarlinSite.abbr
    ]
    headers, row_stub = make_output_stubs(methods)
    methods = sorted(set(methods))
    with args.instl as instl:
        raw_sites = instl.read().split()
    sites, unwrapped = wrap_sites(raw_sites, methods)
    structs = _sites.get_structs([_s[0] for _s in sites], methods)
    with args.outsv as outsv:
        outsv.write(headers)
        for sid, seq in zip(sids, seqs):
            counts = _Counts(seq, structs)
            rows = cbcalc(sid, row_stub, sites, counts, methods)
            outsv.writelines(rows)


if __name__ == "__main__":
    _sys.exit(main())
