"""Written by Kyle Friend at Washington and Lee University in May 2019"""

import re
import copy
import os
import pickle

"""Input is a series of sam alignment files which are used to generate a dictionary
that contains sequence information as well as information about the number of C
to U conversions that were detected during sequencing.
"""
def main():
    Result_Values = {}
    Master_Dict = {}
    (RefSeq_Dict, TempScore_Dict) = FastaReader()
    all_genes = list(RefSeq_Dict.keys())
    all_files = list(os.listdir(r'G:/DecayKinetics'))
    for file in all_files:
        print(file)
        full_path = (r'G:/DecayKinetics/' + file)
        (Master_Dict, TempScore_Dict) = SamReader(RefSeq_Dict, TempScore_Dict, Master_Dict, full_path, file)
    return(Master_Dict)

"""Read input fasta file with transcriptome and generates total mRNA length from
nucleotide sequence. Information is stored in a dictionary"""
def FastaReader():
    Sequences = {}
    Counts = {}
    genes = []
    seq = ''
    seqs = []
    with open(r'E:/DecayKinetics/All_mRNAs.fasta', 'r') as infile:
        for line in infile:
            line = line.rstrip()
            if '>' in line:
                genes.append(line[1:])
                seq = seq.upper()
                seqs.append(seq)
                seq = ''
            else:
                seq = (seq + line)
        seqs.append(seq)
        val = seqs.pop(0)

    for i in range(len(genes)):
        Sequences[genes[i]] = seqs[i]
        Counts[genes[i]] = [0, 0]
    return (Sequences, Counts)

"""Scans through sam alignment files for C to T conversions that occur during sequencing.
The number of reads with a sequence conversion are then normalized to the total number
of reads as a percentage."""
def SamReader(Ref, st_Score, Totals, file_name, sam_ID):
    with open(file_name, 'r') as infile:
        for line in infile:
            line = line.split()
            """line[2] is the gene, and line[3] is the start position. line[9]
            is the sequence for comparison, and line[5] is the CIGAR string."""
            if len(line[9]) > 100 and line[5] != '*':
                seq = Ref[line[2]]
                counts = st_Score[line[2]]
                counts[1] += 1
                if 'I' in line[5] or 'D' in line[5]:
                    start = int(line[3]) - 1
                    end = start + len(line[9])
                    seq = seq[start:end]
                    q_seq = line[9]
                    cig = re.split('I|D|M', line[5])
                    cig.pop(-1)
                    cig = [int(x) for x in cig]
                    high = max(cig)
                    high = cig.index(high)
                    if high == 0:
                        q_seq = q_seq[:cig[0]]
                    elif high == len(cig) - 1:
                        q_seq = q_seq[-cig[-1]:]
                    else:
                        start = sum(cig[:high])
                        end = sum(cig[:high + 1])
                        q_seq = q_seq[start:end]
                    temp = []
                    try:
                        pos = seq.index(q_seq[0:10])
                    except ValueError:
                        continue
                    seq = seq[pos:]
                    for i in range(len(q_seq)):
                        try:
                            if q_seq[i] != seq[i]:
                                temp.append(seq[i] + q_seq[i])
                        except IndexError:
                            break
                else:
                    start = int(line[3]) - 1
                    end = start + len(line[9])
                    seq = seq[start:end]
                    temp = []
                    if line[9] != seq:
                        q_seq = line[9]
                        for i in range(len(q_seq)):
                            if q_seq[i] != seq[i]:
                                temp.append(seq[i] + q_seq[i])
                if len(temp) < 20 and 'TC' in temp:
                    counts[0] += 1
                st_Score[line[2]] = counts
            else:
                continue
        sub_genes = list(st_Score.keys())
        for j in sub_genes:
            counts = st_Score[j]
            if counts[1] > 0:
                counts.append(counts[0]/counts[1]*100)
            else:
                counts.append(0)
            counts = tuple(counts)
            Totals[(j, sam_ID)] = counts
            st_Score[j] = [0, 0]
    infile.close()
    return(Totals, st_Score)

if __name__ == '__main__':
    All_Results = main()
    pickle.dump(All_Results, open(r'E:/DecayKinetics/Kinet_Dict.p', 'wb'))
    #Due to high compute time, dictionary was temporarily stored locally using pickle.
