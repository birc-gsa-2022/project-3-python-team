from __future__ import annotations
import argparse
from typing import Iterable
from radix_sort import *
from typing import NewType
from parsers import parse_fasta, parse_fastq

saidx = NewType("saidx", int)


def main():
    argparser = argparse.ArgumentParser(
        description="Exact matching using a suffix array")
    argparser.add_argument("genome", type=argparse.FileType('r'))
    argparser.add_argument("reads", type=argparse.FileType('r'))
    args = argparser.parse_args()
    # print(f"Find every reads in {args.reads.name} " +
    #      f"in genome {args.genome.name}")

    genome = parse_fasta(args.genome)
    reads = parse_fastq(args.reads)

    out = []
    for chr in genome:
        sa = sa_construction_nsq(genome[chr])
        for read in reads:
            hits = pattern_match(genome[chr], reads[read], sa)
            for hit in hits:
                out.append(
                    f'{read}\t{chr}\t{hit+1}\t{len(reads[read])}M\t{reads[read]}')

    out.sort()
    print('\n'.join(out))


def sa_construction_nsq(x: str) -> list[int]:
    '''
    Suffix array construction in O(n^2), using radix sort


    Put input string into a set O(n) -> sort set to make alphabet string O(n * alphabet)
     -alphabet size is bounded by n anyway.
    Buckets constructed as dict, in order of alphabet string. O(n)
    Input string is split into suffixes, adding $ so each has length n, O(n*n)
     -could get rid of slicing and just use index, but am lazy, and already implemented radix using strings (and we're n^2 anyway)
    radix sort suffixes O(n*n)
    one last O(n) running through sorted suffixes to get their starting index (because lazyness.)
    '''
    n = len(x)
    if n == 0:
        return []
    alpha = alphabet_from_input_string(x)
    suffixes = [x[i:]+i*'$' for i, _ in enumerate(x)]

    sort = radix_sort(suffixes, alpha)
    return [n-len(x.split('$')[0]) for x in sort]


def sa_construction_simple_nsq(x: str) -> list[int]:
    '''
    simple sa construction uing builtin strcmp and sort, for comparisson
    '''
    suf = [x[i:] for i, _ in enumerate(x)]
    suf = sorted(suf)
    return [len(x) - len(i) for i in suf]


def lower_bound(a: str, i: int, lo: int, hi: int, x: str, sa: list[int]) -> int:
    ''' binary search suffix array for lower bound for occurence of character a at position i in a pattern, in string x, given Suffix array SA of x'''
    while hi > lo:
        mid = (hi+lo)//2
        if x[(sa[mid] + i) % len(x)] < a:
            lo = mid+1
        else:
            hi = mid
    return lo


def upper_bound(a: str, i: int, lo: int, hi: int, x: str, sa: list[int]) -> int:
    ''' binary search suffix array for upper bound for occurence of character a at position i in a pattern, in string x, given Suffix array SA of x'''
    while hi > lo:
        mid = (hi+lo)//2
        if ord(x[(sa[mid] + i) % len(x)]) < ord(a)+1:
            lo = mid+1
        else:
            hi = mid
    return lo


def pattern_match(x: str, p: str, sa: list[int]) -> Iterable[int]:
    ''' report all occurences of pattern p in string x, given Suffix array: SA of x. '''
    high = len(sa)
    low = 0

    if len(p) == 0:
        return []
    for i, a in enumerate(p):
        low = lower_bound(a, i, low, high, x, sa)
        high = upper_bound(a, i, low, high, x, sa)

    return sa[low:high]


# Work in progress below
def radix_by_index(x: str, u: list[int], alphabet: list[str]) -> dict[str, str]:
    buckets = {i: 0 for i in alphabet}
    triplets = get_triplet_map(x, u)

    arr = [i for i in u]  # make copy so we don't override u
    sa = [0] * len(u)
    for i in range(3):
        for trip in arr:
            buckets[triplets[trip][i]] += 1
        accsum = 0
        for bucket in buckets:
            buckets[bucket], accsum = accsum, accsum + buckets[bucket]
        for trip in arr:
            sa[buckets[triplets[trip][i]]] = trip
            buckets[triplets[trip][i]] += 1
        for bucket in buckets:
            buckets[bucket] = 0
        arr, sa = sa, arr

    i = 0
    sigma: dict[str, str] = {}
    for trip in arr:
        if triplets[trip] not in sigma:
            sigma[triplets[trip]] = str(i)
            i += 1
    sigma = {str(sigma[key]): key for key in sigma}
    return sigma


def get_triplet_map(x: str, u: list[int]) -> dict[int, str]:
    triplets: dict[int, str] = {}
    for suf in u:
        suffix = x[suf:]
        if len(suffix) >= 3:
            triplet = suffix[:3]
        else:
            triplet = suffix + (3-len(suffix)) * suffix[-1]
        triplets[suf] = triplet
    return triplets


def skew(x: str) -> list[int]:
    sa_0: list[int] = []
    sa_1: list[int] = []
    sa_2: list[int] = []
    sa_12: list[int] = []
    for i, _ in enumerate(x):
        if i % 3 == 0:
            sa_0.append(i)
        elif i % 3 == 1:
            sa_1.append(i)
        else:
            sa_2.append(i)
    sa_12.extend(sa_1)
    sa_12.extend(sa_2)
    alpha = sorted(set(x))

    sigma = radix_by_index(x, sa_12, alpha)
    inverse = {sigma[key]: key for key in sigma}
    if len(sigma) != len(sa_12):
        print(sigma)
        triplets = get_triplet_map(x, sa_12)
        print(triplets)
        sa_12 = skew(''.join([inverse[triplets[trip]] for trip in sa_12]))
    else:
        special_case: list[int] = []
        buckets: dict[str, int] = {}
        print(sigma)
        for i in sigma:
            buckets[i] = 0
        for trip in sa_0:
            if trip + 1 < len(x):
                buckets[x[trip+1]] += 1
            else:
                special_case.append(trip)
        accsum = 0
        for bucket in buckets:
            buckets[bucket], accsum = accsum, accsum + buckets[bucket]
        arr: list[int] = [0] * (len(sa_0)-len(special_case))
        for trip in sa_0:
            if trip + 1 < len(x):
                arr[buckets[x[trip+1]]] = trip
                buckets[x[trip+1]] += 1
        sa_0 = special_case
        sa_0.extend(arr)
    print(sa_0, sa_12)
    # merge
    i, j = 0, 0
    while i < len(sa_12) and j < len(sa_0):
        pass

    return None


if __name__ == '__main__':
    main()
