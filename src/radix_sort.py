
def radix_sort(inp: "list[str]", alphabet="$abcdefghijklmnopqrstuvwxyz") -> "list[str]":
    '''
    this implementation assumes all input strings are of length m.
    also needs an (ordered) alphabet to run, english alphabet + sentinel added as default...
    '''
    m = [i-1 for i in range(len(inp[0]), 0, -1)]
    arr = [None for _ in inp]
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


def alphabet_from_input_string(x: str) -> str:
    s = set(x)
    order = sorted(s)
    return '$' + ''.join(order)
