"""Site wrappers implementing contrast calculation methods.

Contains site wrappers, function for site digitizing, and dumper-loader
functions for wrapped site collections.
"""

import cPickle
import decimal

__all__ = ["digitize", "dump_sl", "load_sl", "get_structs", "Site",
           "MarkovSite", "PevznerSite", "KarlinSite"]

nucls = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

dnucls = {
        'N': [0, 1, 2, 3],
        'A': [0], 'B': [1, 2, 3],
        'C': [1], 'D': [0, 2, 3],
        'G': [2], 'H': [0, 1, 3],
        'T': [3], 'V': [0, 1, 2],
        'M': [0, 1], 'K': [2, 3],
        'R': [0, 2], 'Y': [1, 3],
        'W': [0, 3], 'S': [1, 2],
        }

def _struct(L, U, M):
    if L <= 0:
        return (0, 0, 0)
    if U <= 0 or L <= U:
        return (L, 0, 0)
    return (L, U, M)

def digitize(site, U, M):
    """Return integer representation of a site.
    
    Returns site structure and digitized site, which is a list of integers,
    each one represents one non-degenerate variant of the given site.
    """
    sites = [0]
    begin = 0
    while site[begin] == 'N':
        begin += 1
    end = len(site)
    while site[end - 1] == 'N':
        end -= 1
    for nucl in site[begin:end]:
        if nucl in nucls:
            digit = nucls[nucl]
            sites = [digit + (dsite << 2) for dsite in sites]
        else:
            sites = [digit + (dsite << 2) for dsite in sites
                     for digit in dnucls[nucl]]
    L = end - begin
    U -= begin
    return _struct(L, U, M), sites

def _addN(start, arr_site, init_site, eL, this_step, next_step, U, M):
    if eL == 1:
        return
    for i in range(start, len(arr_site)):
        if arr_site[i] == 'N':
            continue
        arr_site[i] = 'N'
        this_step.append(digitize(arr_site, U, M))
        _addN(i+1, arr_site, init_site, eL-1, next_step, this_step, U, M)
        arr_site[i] = init_site[i]

class Site():
    
    """A parent site wrapper class.
    
    Instance methods:
        get_length -- return int length of the site either considering
                      or ignoring N positions.
        calc_observed -- calculate observed number of site occurences.
        calc_expected -- calculate expected number of site occurences.
        calc_contrast -- calculate ratio of observed to expected number.
    
    Instance variables:
        str_init -- string representation of the site.
        T -- length of str_init.
        str_site -- string representation of the site with gap trimmed
                    in case of bipartite site.
        L -- length of str_site.
        eL -- effective length of the site (without N positions).
        U -- gap position in case of bipartite site, else 0.
        M -- gap length in case of bipartite site, else 0.
        dsite -- digitized site, contains site structure and its digital
                 form.
        struct -- site structure, tuple(L, U, M).
        structs -- a set of all structures required for calculation of
                   contrast and expected number of the site.
    """
    
    def __init__(self, site, maxlen=10):
        """Site(site)
        
            site -- string, could contain ATGCBDHVSWMKRYN symbols.
            maxlen -- site length cutoff, default 10; raise ValueError
                      if site length is greater than the cutoff.
        """
        self.str_init = site.upper().strip('N')
        self.T = len(self.str_init)
        if self.T > 8 and "NNN" in self.str_init:
            uhalf = self.str_init[:6].strip('N')
            dhalf = self.str_init[-6:].strip('N')
            self.str_site = uhalf + dhalf
            self.U = len(uhalf)
        else:
            self.str_site = self.str_init
            self.U = 0
        self.L = len(self.str_site)
        if self.L > maxlen:
            raise ValueError("Site is too long: %d", self.L)
        self.M = self.T - self.L
        self.eL = len(site.upper().replace('N', ''))
        self.dsite = digitize(self.str_site, self.U, self.M)
        self.struct = self.dsite[0]
        self.structs = {self.struct}
        self._prepare()
    
    def _prepare(self):
        """Site preparation for expected number calculation."""
        pass
    
    def get_length(self, effective=True):
        """Get length of the site.
        
            effective -- True (default) if effective length should be
        returned, else False. Effective length is a number of meaning
        (not N) postitions.
        """
        return self.eL if effective else self.L
    
    def calc_observed(self, counts):
        """Calculate observed number of occurences of the site by counts.
        
            counts -- Counts object.
        """
        return counts.get_count(self.dsite)
    
    def calc_expected(self, counts):
        """Estimate expected number of occurences of the site.
        
        Simplest estimation is used based on the assumption that all
        positions of the site are independent.
            counts -- Counts object.
        """
        expected = counts.get_total(self.struct)
        for nucl in self.str_site:
            expected *= counts.get_freq(((1, 0, 0), dnucls[nucl]))
        return expected
    
    def calc_contrast(self, counts):
        """Calculate contrast ratio.
        
            counts -- Counts object.
        """
        expected = self.calc_expected(counts) or float("NaN")
        return self.calc_observed(counts) / expected
    
    def __len__(self):
        return self.T
    
    def __str__(self):
        return self.str_init
    
    def __repr__(self):
        return "<%s '%s' of %d-%d-%d structure>" % (
                self.__class__.__name__, self.str_init,
                self.L, self.U, self.M
                )

class MarkovSite(Site):
    
    """Site wrapper implementing Mmax based expected number calculation.
    
    Override following Site methods:
        calc_expected -- simple Mmax based method implementation (see
                         Schbath et al., 1995; Gelfand et al., 1997).
    """
    
    def _prepare(self):
        self.rpart = digitize(self.str_site[1:], self.U-1, self.M)
        self.structs.add(self.rpart[0])
        self.lpart = digitize(self.str_site[:-1], self.U, self.M)
        self.structs.add(self.lpart[0])
        self.cpart = digitize(self.str_site[1:-1], self.U-1, self.M)
        self.structs.add(self.cpart[0])
    
    def calc_expected(self, counts):
        """Estimate expected number of the site with Mmax based method."""
        if self.eL == 1:
            return counts.get_total((1, 0, 0)) * len(self.dsite[1]) / 4.0
        div = counts.get_count(self.cpart)
        if div == 0:
            return float('NaN')
        num = counts.get_count(self.lpart) * counts.get_count(self.rpart)
        return float(num) / div

class PevznerSite(Site):
    
    """Site wrapper implementing Pevzner's expected number calculation.
    
    Override following Site methods:
        calc_expected -- Pevzner's method implementation (see
                         Pevzner et al., 1989).
    """
    
    def _prepare(self):
        arr_site = list(self.str_site)
        self.singleN = []
        self.doubleN = []
        for i in range(self.L):
            if arr_site[i] == 'N':
                continue
            arr_site[i] = 'N'
            dsite = digitize(arr_site, self.U, self.M)
            self.structs.add(dsite[0])
            self.singleN.append(dsite)
            for j in range(i+1, self.L):
                if arr_site[j] == 'N':
                    continue
                arr_site[j] = 'N'
                dsite = digitize(arr_site, self.U, self.M)
                self.structs.add(dsite[0])
                self.doubleN.append(dsite)
                arr_site[j] = self.str_site[j]
            arr_site[i] = self.str_site[i]
    
    def calc_expected(self, counts):
        """Estimate expected number of the site with Pevzner's method."""
        if self.eL == 1:
            return counts.get_total((1, 0, 0)) * len(self.dsite[1]) / 4.0
        div = 1.0
        for dsite in self.doubleN:
            div *= counts.get_count(dsite)
        if div == 0.0:
            return float('NaN')
        div = pow(div, 2.0 / (self.eL ** 2 - self.eL))
        num = 1.0
        for dsite in self.singleN:
            num *= counts.get_count(dsite)
        return pow(num, 2.0 / self.eL) / div

class KarlinSite(Site):
    
    """Site wrapper implementing Karlin's expected number calculation.
    
    Override following Site methods:
        calc_expected -- Karlin's method implementation (see
                         Karlin et al., 1994).
    """
    def _prepare(self):
        arr_site = list(self.str_site)
        self.oddN = []
        self.evenN = []
        _addN(0, arr_site, self.str_site, self.eL, self.oddN, self.evenN,
              self.U, self.M)
        for struct, dsite in self.oddN:
            self.structs.add(struct)
        for struct, dsite in self.evenN:
            self.structs.add(struct)
    
    def calc_expected(self, counts):
        """Estimate expected number of the site with Karlin's method."""
        if self.eL == 1:
            return counts.get_total((1, 0, 0)) * len(self.dsite[1]) / 4.0
        div = decimal.Decimal(1)
        for dsite in self.evenN:
            div *= decimal.Decimal(counts.get_freq(dsite))
        if div == decimal.Decimal(0):
            return float('NaN')
        num = decimal.Decimal(counts.get_total(self.struct))
        for dsite in self.oddN:
            num *= decimal.Decimal(counts.get_freq(dsite))
        return float(num / div)

def dump_sl(site_list, ounsl):
    """Dump collection of wrapped sites."""
    with ounsl:
        cPickle.dump(site_list, ounsl, -1)

def load_sl(innsl):
    """Load collection of wrapped sites from the dump file."""
    with innsl:
        return cPickle.load(innsl)

def get_structs(site_list):
    """Return a set of structures of all sites and all their subsites.
    
    Use this function to obtain structures for word counts calculation.
    """
    structs = set()
    for site in site_list:
        structs.update(site.structs)
    max_structs = dict()
    for struct in structs:
        l, u, m = struct
        max_structs[(u, m)] = max(l, max_structs.get((u, m), 0))
    structs = set()
    for (u, m), l in max_structs.items():
        structs.add((l, u, m))
    return structs
