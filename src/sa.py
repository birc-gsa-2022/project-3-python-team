from __future__ import annotations
import argparse
from typing import Iterable
from radix_sort import *
from typing import NewType

saidx = NewType("saidx", int)


def main():
    argparser = argparse.ArgumentParser(
        description="Exact matching using a suffix array")
    argparser.add_argument("genome", type=argparse.FileType('r'))
    argparser.add_argument("reads", type=argparse.FileType('r'))
    args = argparser.parse_args()
    print(f"Find every reads in {args.reads.name} " +
          f"in genome {args.genome.name}")


def sa_construction_nsq(x: str) -> "list[int]":
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
    alpha = alphabet_from_input_string(x)
    suffixes = [x[i:]+i*'$' for i, _ in enumerate(x)]

    sort = radix_sort(suffixes, alpha)
    return [n-len(x.split('$')[0])-1 for x in sort]


def sa_construction_simple_nsq(x: str) -> "list[int]":
    '''
    simple sa construction uing builtin strcmp and sort, for comparisson
    '''
    suf = [x[i:] for i, _ in enumerate(x)]
    suf = sorted(suf)
    return [len(x) - len(i) for i in suf]


def lower_bound(a: str, i: int, lo: int, hi: int, x: str, sa: list[int]) -> int:
    while hi > lo:
        mid = (hi+lo)//2
        if x[(sa[mid] + i) % len(x)] < a:
            lo = mid+1
        else:
            hi = mid
    return lo


def upper_bound(a: str, i: int, lo: int, hi: int, x: str, sa: list[int]) -> int:
    while hi > lo:
        mid = (hi+lo)//2
        if ord(x[(sa[mid] + i) % len(x)]) < ord(a)+1:
            lo = mid+1
        else:
            hi = mid
    return lo


def pattern_match(x: str, p: str, sa: "list[int]") -> "Iterable[int]":
    high = len(sa)
    low = 0

    for i, a in enumerate(p):
        low = lower_bound(a, i, low, high, x, sa)
        high = upper_bound(a, i, low, high, x, sa)

    return sa[low:high]


# Work in progress below
def radix_by_index(x: str, u: list[int]) -> list[int]:
    alphabet = sorted(set(x))
    buckets = {i: 0 for i in alphabet}
    sa = [0] * len(u)
    for i in range(3):
        for trip in u:
            buckets[x[trip + 3 - i]] += 1
        accsum = 0
        for bucket in buckets:
            buckets[bucket], accsum = accsum, accsum + buckets[bucket]
        for trip in u:
            sa[buckets[x[trip+3-i]]] = trip
            buckets[x[trip+3-i]] += 1
        for bucket in buckets:
            buckets[bucket] = 0
        u, sa = sa, u

    return


def skew(x: str) -> list[int]:
    sa_0 = []
    sa_12 = []
    for i, _ in enumerate(x):
        if i % 3 == 0:
            sa_0.append(i)
        else:
            sa_12.append(i)


if __name__ == '__main__':
    # main()
    x = "mississippi$"
    foo = sa_construction_nsq(x)
    print(pattern_match(x, "s", foo))
