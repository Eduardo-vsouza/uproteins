import sys
import os
import subprocess

from Bio import SeqIO
import pandas as pd


class SDInspection(object):
    def __init__(self, args, filetype, folder):
        """
        Looks for a Shine-Dalgarno sequence using the script free_align.pl from free2bind package.
        :param args:
        :param filetype: either genome or transcriptome. Lowercase.
        """
        self.rRNA = self.__get_sequences(args.rrna)[::-1][:13]  # one of the RNA strings must be 3'-5'.
        self.filetype = filetype
        self.results = pd.read_csv(f'{folder}/post_perc/{filetype}_results_02.txt', sep='\t')
        self.coordinates = self.results["Genome Coordinates"].tolist()
        self.entries = self.results["Protein"].tolist()
        self.genome = self.__get_sequences(args.genome)
        self.upstreamSequences = self.__extract_upstream()
        self.pyPath = sys.path[0]
        self.freeAlignPath = f'{self.pyPath}/dependencies/free2bind/free_align.pl'

    @staticmethod
    def __get_sequences(fasta):
        records = SeqIO.parse(fasta, 'fasta')
        seqs = [str(record.seq) for record in records][0]
        return seqs

    def __extract_upstream(self):
        """
        Extracts the 21 nucleotides upstream from the START codon of each ORF present in self.results.
        :return:
        """
        upstream_seqs = []
        for coords, entry in zip(self.coordinates, self.entries):
            splat_coords = coords.split("-")
            if 'forward' in entry:
                start = int(splat_coords[0])
                upstream = self.genome[start-22: start]  # index starts at 0. Need the 21 nucleotides upstream
                upstream_seqs.append(upstream)
            elif 'reverse' in entry:
                start = int(splat_coords[1])
                upstream = self.__complement(self.genome[start: start+22][::-1])
                upstream_seqs.append(upstream)
        return upstream_seqs

    @staticmethod
    def __complement(seq):
        """
        :param seq: sequence to generate the reverse complement for
        :return: complementary sequence
        """
        nucs = {"A": "T", "T": "A", "G": "C", "C": "G"}
        compl = ""
        for nuc in seq:
            compl += nucs[nuc]
        return compl

    def check_free_energy(self):
        """
        Runs free_align.pl to check which part of the sequence has the minimum free energy when binding to the rRNA,
        i.e, has the most stable binding. """
        energies = []
        for upstream_seq in self.upstreamSequences:
            cmd = f'{self.freeAlignPath} -e {upstream_seq} {self.rRNA}'
            free_energy = subprocess.check_output(cmd, shell=True)
            energies.append(free_energy)
        print(self.upstreamSequences)
        print(energies)