#! /usr/bin/python

import argparse as ap
import sys

import cbclib.counts as cnt
import cbclib.sites as st


def load_sites(instl, wrapper, len_cutoff=8):
	sites = set()
	with instl:
		for line in instl:
			try:
				site = wrapper(line.strip("\n\t\r Nn-."))
			except ValueError:
				sys.stderr.write("%s is too long, skipped.\n" % str(site))
				continue
			length = site.L
			if length > len_cutoff:
				sys.stderr.write("%s is too long, skipped.\n" % site)
				continue
			sites.add(site)
	structs = st.get_structs(sites)
	return sites, structs

if __name__ == "__main__":
	parser = ap.ArgumentParser(description="Contrast calculation")
	parser.add_argument(
	        "-f", "--fasta", dest="inseq", metavar="file", required=True,
	        help="Input fasta file, may be gzipped."
	        )
	parser.add_argument(
	        "-s", "--sites", dest="instl", metavar="file",
	        type=ap.FileType('r'), default=sys.stdin,
	        help="Input list of sites, one-per-line."
	        )
	parser.add_argument(
	        "-m", "--method", dest="method", default="karlin",
	        choices=["mmax", "pevzner", "karlin"],
	        help="Method of expected frequency calculation, " +
	        "Karlin's method is default."
	        )
	parser.add_argument(
	        "-o", "--out", dest="outsv", metavar="file",
	        type=ap.FileType('w'), default=sys.stdout,
	        help="Output tabular (.tsv) file, default is stdout"
	        )
	args = parser.parse_args()
	wrapper = {"mmax": st.MarkovSite, "pevzner": st.PevznerSite,
	           "karlin": st.KarlinSite}[args.method]
	sites, structs = load_sites(args.instl, wrapper, 10)
	counts = cnt.calc_all(args.inseq, structs)
	title = ("Site\tObserved number\tExpected number (%s)\t" +
	         "Contrast ratio\tTotal number\n") % args.method.capitalize()
	ouline = "{Site}\t{No:d}\t{Ne:.2f}\t{Ratio:.3f}\t{Total:.0f}\n"
	with args.outsv as outsv:
		outsv.write(title)
		for wrapped in sorted(sites, key=str):
			vals = dict()
			vals["Site"] = wrapped.str_init
			vals["Length"] = counts.get_total(wrapped.struct)
			vals["No"] = wrapped.calc_observed(counts)
			vals["Ne"] = wrapped.calc_expected(counts)
			vals["Ratio"] = vals["No"] / (vals["Ne"] or float("NaN"))
			outsv.write(ouline.format(**vals))

