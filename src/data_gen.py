import random
import argparse
import sys
from typing import Callable, Iterable


def main() -> int:
    argparser = argparse.ArgumentParser(
        description="generate fasta or fast q files, with sequences generated by various methods eg. repetetive, random..\nWhen generating more than one sequence, sizes will increas following length * n"
    )
    argparser.add_argument("filetype", type=str, choices=[
                           "fasta", "fastq"], help='specify type of file to generate')
    argparser.add_argument("seq_type", type=str, choices=[
                           "markov", "random", "repetetive", "identical", "palindrome"], help='specify method to generate from')
    argparser.add_argument("-l", "--length", type=int, default=100, nargs='?',
                           help='length of first sequence, default = 100')
    argparser.add_argument("-n", "--n_sequences", type=int, nargs='?',
                           default=10, help='number of sequences to generate of each length, default = 10')
    argparser.add_argument("-r", "--recurrences", type=int, nargs='?', default=10, help='number of sequences of progressivly longer length')
    argparser.add_argument("-o", "--output_file", type=str,
                           default="stdout", help='specify name of output file, default = STDOUT', nargs='?')

    args = argparser.parse_args()
    outdict = foo(args.n_sequences, args.length, args.recurrences, fundir[args.seq_type])
    match args.output_file:
        case "stdout":
            for key in outdict:
                print(f'{typedir[args.filetype]}{key}\n{outdict[key]}')
        case _:
            with open(f'{args.output_file}.fast{"q" if (args.filetype == "fastq") else "a"}', mode='w') as f:
                for key in outdict:
                    f.write(f'{typedir[args.filetype]}{key}\n{outdict[key]}\n')

    return 0


DNA_alphabet = ['A', 'T', 'C', 'G']
DNA_weights = [35, 35, 15, 15]
typedir = {"fasta": ">", "fastq": "@"}

'''
Generating alphabets to test:
'''

def make_alphabet(size: int) -> list:
    '''
    Makes an alphabet string with a maximum size of 50.
    E.g.:
    >>> make_alphabet(4) -> "abcd"
    '''
    letters = 'abcdefghijklmnopqrstuvxyzABCDEFGHIJKLMNOPQRSTUVXYZ'
    alphabet = letters[0:size-1]
    return list(alphabet)


def make_random_weights(alphabet: Iterable) -> list:
    '''
    Creates list of random weights corresponding to an alphabet string. 
    '''
    weights = []
    for i in len(alphabet):
        weights.append(random.randint(0, 100))
    return weights


'''
Generating strings to test:
'''


def one_letter(length: int, alphabet: Iterable[str]) -> str:
    return "A"*length


def tot_rand(length: int, alphabet: Iterable[str]) -> str:
    s = ""
    for i in range(length):
        s += random.choice(alphabet)
    return s


def prob_rand(length: int, alphabet: Iterable[str], weights: list) -> str:
    s = random.choices(alphabet, weights, k=length)
    return "".join(s)


def repeat_rand(length: int, alphabet: Iterable[str], rep_length: int = 10) -> str:
    rep = ""
    for r in range(rep_length):
        rep += random.choice(alphabet)
    seq = rep*(length//rep_length)+rep[0:(length % rep_length)]
    return seq


def repeat_rand_palindrome(length: int, alphabet: Iterable[str], palindrome_length: int = 20) -> str:
    '''
    Returns a sequence of given length with palindromes of total length of palindrome_length.
    E.g.:
    >>> DNA_alphabet=['A','T','C','G']
    >>> repeat_rand_palindrome(20,6,DNA_alphabet) -> CAGGACCAGGACCAGGACCA
    '''
    rep = ""
    for r in range(palindrome_length//2):
        rep += random.choice(alphabet)
    seq = (rep+rep[::-1])*(length//palindrome_length) + \
        (rep+rep[::-1])[0:(length % palindrome_length)]
    return seq


def DNA_markov(length: int, alphabet: Iterable[str]) -> str:
    DNA_alphabet = ['A', 'T', 'C', 'G']
    s = random.choice(DNA_alphabet)
    letter_chances = {"A": [70, 10, 10, 10], "T": [
        10, 70, 10, 10], "C": [15, 15, 55, 15], "G": [15, 15, 15, 55]}
    for _ in range(length-1):
        s = s+"".join(random.choices(DNA_alphabet, letter_chances[s[-1]], k=1))
    return s


fundir: dict[str, Callable[[int, str], str]] = {
    "markov": DNA_markov, "random": tot_rand, "repetetive": repeat_rand, "identical": one_letter, "palindrome": repeat_rand_palindrome}


'''
Write to fasta:
'''


def foo(times: int, length: int, reccurences: int, f: Callable[[int, str], str]) -> "dict[str, str]":
    out: dict[str, str] = {}
    for r in range(reccurences):
        for i in range(times):
            out[f"{f.__name__}_{length*(r+1)}_{i}"] = f(length*(r+1), "ACGT")
    return out


if __name__ == '__main__':
    main()

sys.exit()


def write_to_fasta_one_letter(length: int, times: int):
    f = open("DNA_test_sequences.fasta", "w")
    f.write("")
    f.close
    f = open("DNA_test_sequences.fasta", "a")
    for i in range(times):
        f.write(">one_letter_{}".format(i)+"\n")
        f.write(one_letter(length, DNA_alphabet))
        f.write("\n")


def write_to_fasta_tot_rand(length: int, times: int):
    f = open("DNA_test_sequences.fasta", "w")
    f.write("")
    f.close
    f = open("DNA_test_sequences.fasta", "a")
    for i in range(times):
        f.write(">tot_rand_{}".format(i)+"\n")
        f.write(tot_rand(length, DNA_alphabet))
        f.write("\n")


def write_to_fasta_prob_rand(length: int, times: int):
    f = open("DNA_test_sequences.fasta", "w")
    f.write("")
    f.close
    f = open("DNA_test_sequences.fasta", "a")
    for i in range(times):
        f.write(">prob_rand_{}".format(i)+"\n")
        f.write(prob_rand(length, DNA_alphabet, DNA_weights))
        f.write("\n")


def write_to_fasta_repeat_rand(length: int, times: int):
    f = open("DNA_test_sequences.fasta", "w")
    f.write("")
    f.close
    f = open("DNA_test_sequences.fasta", "a")
    for i in range(times):
        f.write(">repeat_rand_{}".format(i)+"\n")
        f.write(repeat_rand(length, 6, DNA_alphabet))
        f.write("\n")


def write_to_fasta_DNA_markov(length: int, times: int):
    f = open("DNA_test_sequences.fasta", "w")
    f.write("")
    f.close
    f = open("DNA_test_sequences.fasta", "a")
    for i in range(times):
        f.write(">DNA_markov_{}".format(i)+"\n")
        f.write(DNA_markov(length * (i+1)))
        f.write("\n")


def write_to_fastq_one_letter(length: int, times: int):
    f = open("DNA_patterns.fastq", "w")
    f.write("")
    f.close
    f = open("DNA_patterns.fastq", "a")
    for i in range(times):
        f.write("@read_{}".format(i)+"\n")
        f.write(one_letter(length))
        f.write("\n")


def write_to_fastq_tot_rand(length: int, times: int):
    f = open("DNA_patterns.fastq", "w")
    f.write("")
    f.close
    f = open("DNA_patterns.fastq", "a")
    for i in range(times):
        f.write("@read_{}".format(i)+"\n")
        f.write(tot_rand(length, DNA_alphabet))
        f.write("\n")


def write_to_fastq_prob_rand(length: int, times: int):
    f = open("DNA_patterns.fastq", "w")
    f.write("")
    f.close
    f = open("DNA_patterns.fastq", "a")
    for i in range(times):
        f.write("@read_{}".format(i)+"\n")
        f.write(prob_rand(length, DNA_alphabet, DNA_weights))
        f.write("\n")


def write_to_fastq_repeat_rand(length: int, times: int, rep_length: int = 11):
    '''
    NOTICE rep_length ARGUMENT!!!
    '''
    f = open("DNA_patterns.fastq", "w")
    f.write("")
    f.close
    f = open("DNA_patterns.fastq", "a")
    for i in range(times):
        f.write("@read_{}".format(i)+"\n")
        f.write(repeat_rand(length, rep_length, DNA_weights))
        f.write("\n")


def write_to_fastq_repeat_rand_palindrome(length: int, times: int, rep_length: int = 11):
    '''
    NOTICE rep_length ARGUMENT!!!
    '''
    f = open("DNA_patterns.fastq", "w")
    f.write("")
    f.close
    f = open("DNA_patterns.fastq", "a")
    for i in range(times):
        f.write("@read_{}".format(i)+"\n")
        f.write(repeat_rand_palindrome(length, rep_length, DNA_weights))
        f.write("\n")


def write_to_fastq_DNA_markov(length: int, times: int):
    f = open("DNA_patterns.fastq", "w")
    f.write("")
    f.close
    f = open("DNA_patterns.fastq", "a")
    for i in range(times):
        f.write("@read_{}".format(i)+"\n")
        f.write(DNA_markov(length))
        f.write("\n")


if __name__ == '__main__':
    '''
    We should try different combinations of these, but this is the general usage.
    '''
    main()
    #write_to_fasta_DNA_markov(1000, 20)
    #write_to_fastq_DNA_markov(10, 5)
    #write_to_fasta_one_letter(20000, 1)
    #write_to_fastq_one_letter(10, 1)
