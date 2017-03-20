#! /usr/bin/env python

"""Compositional Bias calculation tool."""

import argparse
import sys
from os.path import basename
from multiprocessing import Process, Value, Lock, Pipe

import cbclib.sites
from cbclib.counts import Counts
from cbclib import __version__


def make_output_stubs(methodset, methods=None, nosids=False):
    """Make headers and row stub for output table.

    Arguments:
        methodset -- sorted set of method abbreviations from cbclib.sites
        methods -- a list of method abbreviations to set the order
        nosids -- bool, should the output contain Sequence ID column

    Returns:
        string -- output table headers
        string -- output table row stub to use with str.format()
    """
    if not methods:
        methods = methodset
    headers = "Sequence ID\tSite\tObserved\t"
    row_stub = "{0}\t{1}\t{2:d}\t"
    total_index = 3
    if nosids:
        headers = "Site\tObserved\t"
        row_stub = "{0}\t{1:d}\t"
        total_index = 2
    start = total_index + 1
    method_dict = dict([(m, start + i*2) for i, m in enumerate(methodset)])
    for method in methods:
        headers += "Expected ({0})\tO/E ratio ({0})\t".format(method)
        index = method_dict[method]
        row_stub += "{{{}:.2f}}\t{{{}:.3f}}\t".format(index, index+1)
    headers += "Total\n"
    row_stub += "{%d:.0f}\n" % total_index
    return headers, row_stub


def wrap_sites(raw_sites, methods, maxlen=10):
    """Wrap raw sites with wrappers implementing specified methods.

    Arguments:
        raw_sites -- a list of sites as strings
        methods -- a list of method abbreviations from cbclib.sites
        maxlen (optional) -- site length cutoff, default 10

    Returns:
        list -- a list of tuples, each contains wrapped versions of the
                site for each method
        list -- a list of site structs to calculate counts for
        dict -- a dict of unwrapped raw sites with reasons why they were
                skipped
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
    sites = [wrapped[0] for wrapped in wrapped_sites]
    structs = cbclib.sites.get_structs(sites, methods)
    return wrapped_sites, structs, unwrapped_sites


def cbcalc(sites, counts):
    """Calculate values for the sites by the counts.

    Arguments:
        sites -- a list of wrapped sites obtained with wrap_sites()
        counts -- cbclib.counts.Counts instance

    Returns:
        list - a list of tuples of calculated values
    """
    vals = []
    for wrapped in sites:
        row_vals = []
        wsite = wrapped[0]
        row_vals.append(wsite.str_init)
        obs = wsite.calc_observed(counts)
        row_vals.append(obs)
        row_vals.append(counts.get_total(wsite.struct_hash))
        for wsite in wrapped:
            exp = wsite.calc_expected(counts)
            row_vals.append(exp)
            row_vals.append(obs / (exp or float("NaN")))
        vals.append(tuple(row_vals))
    return vals


class CBcalcProcess(Process):
    def __init__(self, index, seqs, sites, structs, locked_pipe):
        super(CBcalcProcess, self).__init__()
        self.index = index
        self.seqs = seqs
        self.sites = sites
        self.structs = structs
        self.oulock, self.oupipe = locked_pipe

    def run(self):
        while self.index.value < len(self.seqs):
            with self.index.get_lock():
                index = self.index.value
                self.index.value += 1
            sid, fasta = self.seqs[index]
            counts = Counts(fasta, self.structs)
            vals = cbcalc(self.sites, counts)
            with self.oulock:
                self.oupipe.send((sid, vals))
        with self.oulock:
            self.oupipe.send(("", []))


def main(argv=None):
    """Main function.

    Parse CLI arguments and execute `cbcalc` function.
    """
    parser = argparse.ArgumentParser(
        description="CBcalc - Compositional Bias calculation tool.",
        add_help=False, usage="\n    %(prog)s ".join([
            "", "--help", "--version",
            "[-s LIST] [-o FILE] [-m N] [-BMPK] FASTA [FASTA ...]",
            "[-s LIST] [-o FILE] [-m N] [-BMPK] PATH -i LIST"
        ])
    )
    parser.add_argument(
        "inseq", metavar="FASTA", nargs="+",
        help="Fasta file, may be gzipped."
    )
    parser.add_argument(
        "-i", "--id", dest="insids", metavar="LIST",
        type=argparse.FileType("r"), help="""File with a
        whitespace-separated list of sequence IDs; makes CBcalc to treat
        the first (and the only in the case) positional argument as the
        PATH stub."""
    )
    parser.add_argument(
        "_none", metavar="PATH", nargs="?", #to make up the help message
        help="""Fasta files path stub with '{}' placeholder for
        sequence IDs, which are listed in the file specified with
        -i/--id option."""
    )
    parser.add_argument(
        "-s", "--site", dest="instl", metavar="LIST",
        type=argparse.FileType('r'), default=sys.stdin, help="""File
        with a whitespace-separated list of sites, default is STDIN."""
    )
    parser.add_argument(
        "-o", "--out", dest="outsv", metavar="FILE",
        type=argparse.FileType('w'), default=sys.stdout,
        help="Output tabular (.tsv) file, default is STDOUT."
    )
    parser.add_argument(
        "-m", "--threads", dest="proc_num", metavar="N", type=int,
        default=1, help="""Number of subprocesses to use, should not
        exceed the number of the input fasta files; default is 1."""
    )
    method_group = parser.add_argument_group(
        description="""Folowing arguments determine which methods of
        expected frequency calculation to use. The order of the
        arguments defines the column order of the output file. If
        no method is specified, default set (mmax, pevzner, burge)
        will be used."""
    )
    method_group.add_argument(
        "-B", "--bernoulli", dest="methods", action="append_const",
        const=cbclib.sites.Site.abbr,
        help="Method based on Bernoulli model."
    )
    method_group.add_argument(
        "-M", "--mmax", dest="methods", action="append_const",
        const=cbclib.sites.MarkovSite.abbr,
        help="Method based on maximum order Markov chain model."
    )
    method_group.add_argument(
        "-P", "--pevzner", dest="methods", action="append_const",
        const=cbclib.sites.PevznerSite.abbr,
        help="Method of Pevzner et al."
    )
    method_group.add_argument(
        "-K", "--burge", dest="methods", action="append_const",
        const=cbclib.sites.KarlinSite.abbr,
        help="Method of Burge et al."
    )
    parser.add_argument(
        "-v", "--version", action="version",
        version=("CBcalc v%s" % __version__), help=argparse.SUPPRESS
    )
    parser.add_argument(
        "-h", "-?", "-help", "--help", action="help",
        help=argparse.SUPPRESS
    )
    args = parser.parse_args(argv)
    if args.insids is not None:
        with args.insids as insids:
            sids = sorted(set(insids.read().split()))
        seqs = [args.inseq[0].format(sid) for sid in sids]
    else:
        seqs = args.inseq
        sids = [basename(seq).split(".fasta")[0] for seq in seqs]
    methods_order = args.methods or [
        cbclib.sites.MarkovSite.abbr,
        cbclib.sites.PevznerSite.abbr,
        cbclib.sites.KarlinSite.abbr
    ]
    methodset = sorted(set(methods_order))
    headers, row_stub = make_output_stubs(methodset, methods_order)
    with args.instl as instl:
        raw_sites = instl.read().split()
    sites, structs, _unwrapped = wrap_sites(raw_sites, methodset)
    with args.outsv as outsv:
        outsv.write(headers)
        fastas = zip(sids, seqs)
        if args.proc_num > 1:
            reciver, sender = Pipe(False)
            locked_sender = (Lock(), sender)
            index = Value("i", 0)
            working = args.proc_num
            for _i in range(working):
                proc = CBcalcProcess(
                    index, fastas, sites, structs, locked_sender
                )
                proc.start()
            while working:
                reciver.poll(None)
                sid, vals = reciver.recv()
                if sid:
                    outsv.writelines([
                        row_stub.format(sid, *row) for row in vals
                    ])
                else:
                    working -= 1
            reciver.close()
        else:
            for sid, seq in fastas:
                counts = Counts(seq, structs)
                vals = cbcalc(sites, counts)
                outsv.writelines([
                    row_stub.format(sid, *row) for row in vals
                ])


if __name__ == "__main__":
    sys.exit(main())
