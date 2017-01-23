"""A compositional bias calculation library.

The library contains all you need to estimate the representation of certain
short words (oligonucleotides) in a given nucleotide text, e.g. genome
sequence, with a set of methods. All of them imply calculation of the ratio
of observed frequency of a word to the frequency one could expect to
observe in case the word is not subjected to any selective pressure for or
against. The ratio is called compositional bias here. The methods of
compositional bias calculation differentiate from each other only by the
technique of estimation of the expected frequency.

There is a word (site) wrapper class for each of the methods, which
implements the corresponding technique of estimation of the expected
frequency, in the `sites` module. The parent wrapper class called `Site`
implements the naive Bernoulli model-based approach.

Use of some methods require extensive calculation of observed frequencies
of subpatterns of the word of interest. That is why, there is the
preliminary stage of calculation of the observed numbers (counts) for all
possible words of a given length. All necessary means for the stage are in
counts module.

The pipeline sketch:

[RAW SITES] --(wrapping with Site subclasses)--> [WRAPPED SITES]
                                                             |
          (listing subword structs to calculate counts for)--|
                         |                                   |
[TEXT] --(calculating subword counts)--> Counts obj          |
                                             |               |
             [OUTPUT] <--(calculating compositional biases)--'
"""

__all__ = ["counts", "sites", "version"]

from .version import full_version as __version__
