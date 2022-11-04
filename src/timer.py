import time
import argparse
import numpy as np
from parsers import parse_fasta, parse_fastq
from sa import sa_construction_simple_nsq, pattern_match
import csv

def main():
    argparser = argparse.ArgumentParser(
        description="Exact matching using a suffix array")
    argparser.add_argument("genome", type=argparse.FileType('r'))
    argparser.add_argument("reads", type=argparse.FileType('r'))
    argparser.add_argument("-o", "--output",type=str,default="out")
    args = argparser.parse_args()
    print(f"Find every reads in {args.reads.name} " +
          f"in genome {args.genome.name}")

    genome = parse_fasta(args.genome)
    reads = parse_fastq(args.reads)
    times=time_construction_and_search(genome, reads)
    with open("{}genome.csv".format(args.output),"w") as f:
        write=csv.writer(f)
        write.writerow(times[0])
    with open("{}reads.csv".format(args.output),"w") as f:
        write=csv.writer(f)
        write.writerows(times[1])

def time_construction_and_search(genome, reads):
    out1=[]
    out2=[] 
    for chr in genome:
        print(len(genome[chr]))
        temp=[]
        for i in range(1000):
            t=time.process_time()
            sa = sa_construction_simple_nsq(genome[chr])
            temp.append(time.process_time()-t)
        out1.append(np.sum(temp))
        read_times=[]
        for read in reads:
            temp=[]
            for i in range(1000):
                t=time.process_time()
                hits = pattern_match(genome[chr], reads[read], sa)
                temp.append(time.process_time()-t)
            read_times.append(np.sum(temp))
        out2.append(read_times)
    return out1, out2

if __name__ == '__main__':
    main()
