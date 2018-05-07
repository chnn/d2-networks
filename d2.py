from functools import reduce
from itertools import product
from math import sqrt


def d2(seq_a, seq_b, k):
    """Compute the d2 statistic for two genomic sequences.

    Parameters
    ----------
    seq_a : str
    seq_b : str
    k : int

    Returns
    -------
    float

    """
    a_counts = count_kmers(seq_a, k)
    b_counts = count_kmers(seq_b, k)

    return reduce(
        lambda a, kmer: a + (a_counts.get(kmer, 0) * b_counts.get(kmer, 0)),
        kmers(seq_a, k),
        0
    )


def d2star(seq_a, seq_b, k):
    """Compute the $d_2^*$ statistic for two genomic sequences.

    This implementation is simple, but far too slow to be of practical use.

    Parameters
    ----------
    seq_a : str
    seq_b : str
    k : int

    Returns
    -------
    float

    """
    ps = letter_probabilities(seq_a + seq_b)
    n_bar = len(seq_a) - k + 1
    m_bar = len(seq_b) - k + 1
    seq_a_counts = count_kmers(seq_a, k)
    seq_b_counts = count_kmers(seq_b, k)
    result = 0

    for w in all_kmers(k):
        p_w_hat = reduce(lambda a, x: a * ps[x.lower()], w, 1)
        x_w_tilde = seq_a_counts.get(w, 0) - (n_bar * p_w_hat)
        y_w_tilde = seq_b_counts.get(w, 0) - (m_bar * p_w_hat)
        result += (x_w_tilde * y_w_tilde) / (sqrt(n_bar * m_bar) * p_w_hat)

    return result


def letter_probabilities(seq):
    """Estimate the probability that a character appears in `seq`.

    Returns
    -------
    dict
        With a key for each character (case-insensitive) appearing in `seq` and
        values being sample probability that this character appears in `seq`.
    """
    ps = {}

    for char in seq:
        ps[char.lower()] = ps.get(char.lower(), 0) + 1

    return {char: count / len(seq) for char, count in ps.items()}


def count_kmers(seq, k):
    """Get counts of all kmers in `seq`.

    Returns
    -------
    dict
        Keys are (lowercased) kmers that appear in `seq`, and values are
        integers.

    """
    counts = {}

    for kmer in kmers(seq, k):
        counts[kmer] = counts.get(kmer, 0) + 1

    return counts


def kmers(seq, k):
    """Get an iterable for all (lowercased) kmers in `seq`.

    """
    for i in range(len(seq) - k + 1):
        yield seq[i:i + k].lower()


def all_kmers(k):
    """Get an iterable for all possible kmers generated from the alphabet
    `['a', 'c', 'g', 't']`.

    """
    return map(lambda x: "".join(x), product('acgt', repeat=k))
