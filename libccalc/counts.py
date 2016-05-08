"""A wrapper and handlers for word counts.

Contains a wrapper for word counts calculated with count_words, which is
either stand-alone util or C dynamic library function accessed through
ctypes.
"""

import os
import ctypes

__all__ = ["Counts", "calc_all"]

class Counts():
	
	"""Wrapper for word counts."""
	
	def __init__(self):
		"""Counts(length, counts, totals)
		
		Counts wrapper constructor. Creates empty Counts object.
		Note: Use calc() and load() to calculate or upload counts.
		"""
		self.counts = dict()
		self.ptrs = []
	
	def free_counts(self):
		if self.ptrs:
			module_path = os.path.dirname(__file__) or "."
			cw_lib = ctypes.CDLL(module_path + "/dll/count_words.so")
			cw_lib.free_counts.argtypes = [ctypes.c_void_p]
			cw_lib.free_counts.restype = None
			for ptr in self.ptrs:
				cw_lib.free_counts(ptr)
			self.ptrs = []
			self.counts.clear()
	
	def get_count(self, site):
		"""Get count of given site.
		
		  site -- (struct, dsites) tuple.
		    struct -- (len, pos, gap) tuple.
		      len -- length of the site.
		      pos -- length of upstream part in case of bipartite site,
		             0 in case of short site.
		      gap -- length of gap in case of bipartite site,
		             0 in case of short site.
		    dsites -- list of digitized sites (see sites module for
		              details).
		
		Return sum of dsite counts.
		"""
		count = 0
		struct, dsites = site
		counts, total = self.counts[struct]
		for dsite in dsites:
			count += counts[dsite]
		return count
	
	def get_nucl_count(self, nucl):
		count = 0
		counts, total = self.counts[(1, 0, 0)]
		for dnucl in nucl:
			count += counts[dnucl]
		return count
	
	def get_freq(self, site):
		"""Get frequency of given site.
		
		  site -- (struct, dsites) tuple.
		    struct -- (len, pos, gap) tuple.
		      len -- length of the site.
		      pos -- length of upstream part in case of bipartite site,
		             0 in case of short site.
		      gap -- length of gap in case of bipartite site,
		             0 in case of short site.
		    dsites -- list of digitized sites (see sites module for
		              details).
		
		Return sum of dsite counts normalized with total number of sites.
		"""
		count = 0.0
		struct, dsites = site
		counts, total = self.counts[struct]
		for dsite in dsites:
			count += counts[dsite]
		return count / total
	
	def get_nucl_freq(self, nucl):
		count = 0.0
		counts, total = self.counts[(1, 0, 0)]
		for dnucl in nucl:
			count += counts[dnucl]
		return count / total
	
	def get_total(self, struct):
		return self.counts[struct][1]
	
	def load(self, incnt):
		"""Upload counts from file.
		
		Upload counts from file obtained with count_words util.
		"""
		with incnt:
			counts = None
			for line in incnt:
				if line.startswith("#"):
					if counts:
						self.counts[struct] = (counts, total)
					vals = line.strip("\n\r#").split(", ")
					title = dict(map(lambda x: x.split("="), vals))
					total = int(title["total_number"])
					site_len = int(title["word_length"])
					half_len = int(title.get("gap_position", 0))
					gap_len = int(title.get("gap_length", 0))
					struct = (site_len, half_len, gap_len)
					num = 1 << site_len * 2
					CountsArr = ctypes.c_long * num
					counts = CountsArr()
					site = 0
				else:
					counts[site] = int(line.split('\t')[-1])
					site += 1
			self.counts[struct] = (counts, total)
	
	def calc(self, fasta_path, struct):
		"""Calculate counts for sequence in the fasta file.
		
		Calculate counts for all words of the length and all shorter ones.
		  fasta_path -- fasta file, may be gzipped.
		  struct -- structure of the words to calculate counts for.
		
		Note: requires count_words/count_words.so library.
		"""
		slen, hlen, glen = struct
		num = slen - hlen
		restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_long) * num)
		module_path = os.path.dirname(__file__) or "."
		cw_lib = ctypes.CDLL(module_path + "/dll/count_words.so")
		if hlen:
			count_words = cw_lib.count_bipart_words
			count_words.restype = restype
			counts = count_words(fasta_path, slen, hlen, glen).contents
		else:
			count_words = cw_lib.count_short_words
			count_words.restype = restype
			counts = count_words(fasta_path, slen).contents
		num_index = 1 << slen * 2
		struct_ = [slen, hlen, glen]
		for i in range(num):
			self.counts[tuple(struct_)] = (counts[i], counts[i][num_index])
			self.ptrs.append(counts[i])
			num_index >>= 2
			struct_[0] -= 1
	
	def __enter__(self):
		return self
	
	def __exit__(self, exc_type, exc_value, traceback):
		self.free_counts()

def calc_all(fasta_path, structs):
	counts = Counts()
	for struct in structs:
		counts.calc(fasta_path, struct)
	return counts
