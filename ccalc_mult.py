#! /usr/bin/python

import argparse as ap
import sys

import libccalc.counts as cnt
import libccalc.sites as st


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
	        "-a", "--acv", dest="inacv", metavar="file", required=True,
	        type=ap.FileType('r'), help="input file with ACv list (.acv)"
	        )
	parser.add_argument(
	        "-i", "--inpath", dest="inpath", metavar="path", required=True,
	        help=".fasta[.gz] files path with %%s placeholder that will "
	        "be replaced with ACv from the input ACv list, only parent " +
	        "folder could be specified, '/%%s.fasta.gz' tag will be " +
	        "added in the case"
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
	        "-o", "--oupath", dest="oupath", metavar="path", required=True,
	        help="output .tsv files path with %%s placeholder or " +
	        "output folder ('/%%s.tsv' tag will be added)"
	        )
	args = parser.parse_args()
	wrapper = {"mmax": st.MarkovSite, "pevzner": st.PevznerSite,
	           "karlin": st.KarlinSite}[args.method]
	sites, structs = load_sites(args.instl, wrapper, 12)
	with args.inacv as inacv:
		acvs = inacv.read().strip().split('\n')
	inpath = args.inpath
	if "%s" not in inpath:
		inpath += "/%s.fasta.gz"
	oupath = args.oupath
	if "%s" not in oupath:
		oupath += "/%s.tsv"
	title = ("Site\tObserved number\tExpected number (%s)\t" +
	         "Contrast ratio\tTotal number\n") % args.method.capitalize()
	ouline = "{Site}\t{No:d}\t{Ne:.2f}\t{Ratio:.3f}\t{Total:.0f}\n"
	for acv in acvs:
		counts = cnt.calc_all(inpath % acv, structs)
		with open(oupath % acv, 'w') as outsv:
			outsv.write(title)
			for wrapped in sorted(sites, key=str):
				vals = dict()
				vals["Site"] = wrapped.str_init
				vals["Total"] = counts.get_total(wrapped.struct)
				vals["No"] = wrapped.calc_observed(counts)
				vals["Ne"] = wrapped.calc_expected(counts)
				vals["Ratio"] = vals["No"] / (vals["Ne"] or float("NaN"))
				outsv.write(ouline.format(**vals))

