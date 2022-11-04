from dataclasses import dataclass


def radix_sort(inp: "list[suffix]", alphabet="$abcdefghijklmnopqrstuvwxyz") -> "list[suffix]":
    '''
    this implementation assumes all input strings are of length m.
    also needs an (ordered) alphabet to run, english alphabet + sentinel added as default...
    '''
    m = [i-1 for i in range(len(inp[0]), 0, -1)]
    arr: list[None | suffix] = [None for _ in inp]
    buckets = {i: 0 for i in alphabet}

    for j in m:
        for bucket in buckets:
            buckets[bucket] = 0
        for x in inp:
            buckets[x[j]] += 1
        accsum = 0
        for bucket in buckets:
            buckets[bucket], accsum = accsum, accsum + buckets[bucket]
        for x in inp:
            arr[buckets[x[j]]] = x
            buckets[x[j]] += 1
        arr, inp = inp, arr
    return inp


@dataclass
class suffix:
    x: str
    i: int

    def __getitem__(self, j: int) -> str:
        if not self.x or self.i + j >= len(self.x):
            return '$'
        return self.x[(self.i+j) % len(self.x)]

    def __len__(self) -> int:
        return len(self.x)

    def __str__(self) -> str:
        return self.x[self.i:]


def radix_bedre(inp: list[int], x: str) -> list[int]:
    pass


def alphabet_from_input_string(x: str) -> str:
    s = set(x)
    order = sorted(s)
    return '$' + ''.join(order)
