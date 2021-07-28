import sys
import os
import subprocess

from Bio import SeqIO
import pandas as pd

from ..sequtils.orflib import ORF, ORFCollection
from ..sequtils.locus import StringTieGFF
from ..sequtils.transcriptomics import TranscriptExtractor


class SDInspection(object):
    def __init__(self, args, filetype, folder, alternatives, subset="Genome", transcriptome_gff=None, transcripts=None):
        """
        Looks for a Shine-Dalgarno sequence using the script free_align.pl from free2bind package.
        :param args: uProteInS arguments
        :param filetype: either genome or transcriptome. Lowercase.
        :param alternatives: the dictionary containing alternative START codons for a given STOP codon, after extending
        it with AltORF method extend_orfs(), preferably.
        """
        self.subset = subset
        self.tORFs = {}
        self.__get_transcripts(transcriptome_gff, assembly=transcripts)
        self.rRNA = self.__get_sequences(args.rrna)[::-1][:13]  # one of the RNA strings must be 3'-5'.
        self.filetype = filetype
        self.folder = folder
        # self.results = pd.read_csv(f'{folder}/post_perc/{filetype}_results_02.txt', sep='\t')
        # self.coordinates = self.results["Genome Coordinates"].tolist()
        # self.entries = self.results["Protein"].tolist()
        self.genome = self.__get_sequences(args.genome)
        # self.upstreamSequences = self.__extract_upstream()
        self.pyPath = sys.path[0]
        self.freeAlignPath = f'{self.pyPath}/dependencies/free2bind/free_align.pl'

        self.alternatives = alternatives
        self.alternatives = self.__extract_upstream()

    def __get_transcripts(self, transcriptome_gff, assembly):
        """

        :param transcriptome_gff: GFF containing the StringTie assembly
        :param transcripts: fasta file containing the transcripts
        :return:
        """
        if transcriptome_gff is not None and assembly is not None:
            self.subset = "Transcriptome"
            gff = StringTieGFF(gff=transcriptome_gff)
            torfs = gff.get_dict()
            sequences = TranscriptExtractor(assembly=assembly)
            transcripts = sequences.get_transcripts()
                        # self.transcripts = transcripts
            for i in torfs:
                orf = torfs[i]
                orf.transcript = transcripts[orf.name]
                self.tORFs[orf.name] = orf
        return self

    @staticmethod
    def __get_sequences(fasta):
        records = SeqIO.parse(fasta, 'fasta')
        seqs = [str(record.seq) for record in records][0]
        return seqs

    def __extract_upstream_old(self):
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
                # print(coords, upstream)
                upstream_seqs.append(upstream)
        return upstream_seqs

    def __extract_upstream(self):
        """
        Extracts the 21 nucleotides upstream from the START codon of each ORF present in a ORFCollection().
        :return: Dictionary with ORFCollections containing the alternative START codons for a given STOP.
        """
        upstream_seqs = {}
        for stop in self.alternatives:
            # print(stop)
            # print(self.alternatives[stop])
            # print([alt.name for alt in self.alternatives[stop]])
            # print(self.alternatives[stop])
            for alt in self.alternatives[stop]:
                # if alt.name != "Discard":
                    # print(alt.name, alt.strand, alt.start)
                if alt.strand == 'forward':
                    if self.subset == "Genome":
                        upstream = self.genome[alt.start-22: alt.start]
                    else:
                        if alt.start-22 > 0:
                            upstream = alt.transcript[alt.start-22: alt.start]
                        else:
                            upstream = ''
                    alt.upstream = upstream
                elif alt.strand == 'reverse':
                    upstream = self.__complement(self.genome[alt.start: alt.start+22][::-1])
                    alt.upstream = upstream
                if self.subset == 'Genome':
                    if str(alt.end) not in upstream_seqs:
                        upstream_seqs[str(alt.end)] = ORFCollection()
                        upstream_seqs[str(alt.end)].add_orf(alt)
                    elif str(alt.end) in upstream_seqs:
                        # print(upstream_seqs[alt.end])
                        upstream_seqs[str(alt.end)].add_orf(alt)
                elif self.subset == "Transcriptome":
                    identifier = f'{alt.transcriptName}_{alt.end}'
                    if identifier not in upstream_seqs:
                        upstream_seqs[identifier] = ORFCollection()
                        upstream_seqs[identifier].add_orf(alt)
                    elif identifier in upstream_seqs:
                        upstream_seqs[identifier].add_orf(alt)
                # else:
                #     if str(alt.end) in upstream_seqs:
                #         upstream_seqs[str(alt.end)].add_orf(alt)
                #     else:
                #         upstream_seqs[str(alt.end)] = ORFCollection()
                #         upstream_seqs[str(alt.end)].add_orf(alt)
        # for stop in upstream_seqs:
            # for orf in upstream_seqs[stop]:
            #     print('upstream', orf.name, orf.strand, orf.end)
            # print([('upstream', orf.name, orf.start, orf.end) for orf in upstream_seqs[stop]])

        return upstream_seqs

    def __define_identifier(self, orf_end, transcript_name=None):
        if self.subset == "Genome":
            identifier = orf_end
        else:
            identifier = f'{transcript_name}_{orf_end}'
        return identifier

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

    def get_free_energy(self):
        alts = {}
        for stop in self.alternatives:
            # print(stop)
            # print(self.alternatives[stop])

            for alt in self.alternatives[stop]:
                # if alt.name != "Discard":
                print(alt)
                print(alt.upstream)
                if alt.upstream != '':

                    cmd = f'{self.freeAlignPath} -e {alt.upstream} {self.rRNA}'
                    energy = subprocess.check_output(cmd, shell=True).strip().rstrip()
                    # print(energy)
                    alt.freeEnergy = float(energy)
                    sd_seq = self.__check_rbs(energy)
                    alt.shineDalgarno = sd_seq
                    # else:
                    #     alt.freeEnergy = 0
                    #     alt.shineDalgarno = "Absent"
                else:
                    alt.freeEnergy = 0
                    alt.shineDalgarno = self.__check_rbs(0)
                identifier = self.__define_identifier(alt.end, transcript_name=alt.transcriptName)
                if identifier not in alts:
                    alts[identifier] = ORFCollection()
                    alts[identifier].add_orf(alt)
                else:
                    alts[identifier].add_orf(alt)
        # for alt in alts:
        #     for orf in alts[alt]:
        #         if orf.name != "Discard":
        #         print(orf.name, orf.start, orf.end)
        #             print(orf.freeEnergy, orf.shineDalgarno)
        return alts

    def get_free_energy_old(self):
        """
        Runs free_align.pl to check which part of the sequence has the minimum free energy when binding to the rRNA,
        i.e, has the most stable binding. """
        energies = []
        sd = []
        for upstream_seq in self.upstreamSequences:
            cmd = f'{self.freeAlignPath} -e {upstream_seq} {self.rRNA}'
            energy = subprocess.check_output(cmd, shell=True).strip().rstrip()
            energies.append(f'{energy}')
            sd_seq = self.__check_rbs(energy)
            sd.append(sd_seq)
        self.results.insert(5, "Free Energy", energies)
        self.results.insert(6, "Shine-Dalgarno", sd)
        return self

    @staticmethod
    def __check_rbs(rbs):
        if float(rbs) >= 0:
            sd_seq = "Leaderless"
        elif -8.4 < float(rbs) < 0:
            sd_seq = "Present. Low to moderate binding."
        elif float(rbs) <= -8.4:
            sd_seq = "Present. Strong binding."
        else:
            sd_seq = "Leaderless"
        return sd_seq

    def save_data(self):
        out = f'{self.folder}/post_perc/{self.filetype}_results_03.txt'
        self.results.to_csv(out, sep='\t', index=False)
        return self


