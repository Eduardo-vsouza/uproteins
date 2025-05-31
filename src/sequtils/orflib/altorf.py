# Copyright © 2021-2025 Eduardo Vieira de Souza
# Copyright © 2021-2025 Adriana Canedo
# Copyright © 2021-2025 Cristiano Valim Bizarro
# Copyright © 2025 Bruno Maestri A Becker
#
# This file is part of uProteInS.
#
# uProteInS is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# uProteInS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# uProteInS. If not, see <https://www.gnu.org/licenses/>.


import pandas as pd
from Bio import SeqIO

from . import ORF, ORFCollection
from ..conversion import Translator
from ..locus import StringTieGFF
from ..transcriptomics import TranscriptExtractor


class AltCodons(object):
    def __init__(self, file, genome, maxsize, subset="Genome", transcriptome_gff=None, assembly=None, testing=False):
        """ I hate this code """
        self.subset = subset
        self.testing = testing
        self.tORFs = self.__check_transcriptome(transcriptome_gff, assembly)
        self.maxsize = maxsize
        self.df = pd.read_csv(file, sep='\t')
        self.__check_testing()
        self.coordinates = self.df["Genome Coordinates"].tolist()
        self.names = self.df["Protein"].tolist()

        self.proteinSequences = self.df["ORF Sequence"].tolist()
        # for i in self.proteinSequences:
            # if i == 'nan':
            # print(i)

        # print(self.proteinSequences)
        self.__genome_records = SeqIO.parse(genome, 'fasta')
        self.genome_seq = [str(record.seq) for record in self.__genome_records]

        self.alternatives = self.__fetch_orfs()
        self.alternatives = self.__extract_peptides()
        # for stop in self.alternatives:
        #     for alt in self.alternatives[stop]:
        #         print(alt.MSPeptides, alt.name)
        #         break

    def __check_testing(self):
        if self.testing:
            self.df = self.df.head(20)

    def __check_transcriptome(self, gff, transcripts_fasta):
        if gff is not None and transcripts_fasta is not None:
            self.subset = "Transcriptome"
            assembly = TranscriptExtractor(assembly=transcripts_fasta)
            transcripts = assembly.get_transcripts()
            # print(transcripts)
            gff = StringTieGFF(gff)
            orf_dict = gff.get_dict()
            # print(orf_dict)
            for gene in orf_dict:
                orf = orf_dict[gene]
                # print(orf, orf.name)

                # orf.transcript = transcripts[orf.name]
                if 'uproteins' in orf.name:
                    rna = '.'.join(orf.name.split(".")[:2])
                    print(rna)
                elif 'gene' in orf.name:
                    rna = orf.name[5:]
                elif 'rna' in orf.name:
                    rna = orf.name[4:]
                else:
                    rna = orf.name
                orf.transcript = transcripts[rna]

            return orf_dict

    def __split_coords(self, i):
        splat = self.coordinates[i].split("-")
        if 'reverse' in self.names[i]:
            start = splat[1]
            end = splat[0]
            strand = 'reverse'
        else:
            start = splat[0]
            end = splat[1]
            strand = 'forward'
        return int(start), int(end), strand

    def __get_transcript_name(self, entry):
        # print(entry)
        gene = entry.split("_")[1]
        if 'gene' in gene:
            # name = gene[5:]
            name = gene
        elif 'rna' in gene:
            # name = gene[4:]
            name = gene
        else:
            splat = gene.split(".")
            # name = f'{splat[0]}.{splat[1]}'
            name = '.'.join(splat[:2])
        return name

    def __get_transcript_coordinates(self, entry):
        splat = entry.split("_")
        coords = splat[len(splat)-2].split("-")
        start, end = coords[0], coords[1]
        # print('transcript coordinates', entry, start, end)
        return int(start), int(end)

    def __fetch_orfs(self):
        """
        :returns a dictionary containing all ORFs with alternative START codons for a given STOP codon.
        """
        alt_check = {}
        alternatives = {}
        print('fetching orfs')
        for i in range(len(self.names)):
            if self.subset == "Genome":
                name = self.names[i]
                start, end, strand = self.__split_coords(i)
                transcript = None
            elif self.subset == "Transcriptome":
                name = self.__get_transcript_name(self.names[i])

                start, end = self.__get_transcript_coordinates(self.names[i])
                # print(name)
                # print(start, end)
                # start, end, strand = self.tORFs[name].start, self.tORFs[name].end, 'forward'
                strand = 'forward'
                # print(name)
                # print(self.tORFs)
                # print(self.tORFs)
                # print(name)
                transcript = self.tORFs[name].transcript
                # print(transcript)
            if strand == 'forward':
                if self.subset == "Genome":
                    seq = self.genome_seq[0][start -1: end]
                else:
                    seq = transcript[start -1: end]
            else:
                seq = self.genome_seq[0][end-1: start][::-1]
                to_comp = Translator(seq)
                seq = to_comp.complement()


            orf = ORF(name=self.names[i], start=int(start), end=int(end), seq=seq, strand=strand, protein_sequence=self.proteinSequences[i])
            orf.transcript = transcript
            orf = self.__fetch_codons(orf)
            orf.transcriptName = name
            identifier = self.__define_identifier(end, transcript_name=orf.transcriptName)
            if identifier not in alt_check:
                alt_check[identifier] = []
            # if start not in alt_check[end]:
            #     alt_check[end].append(start)
            if identifier not in alternatives:
                alternatives[identifier] = ORFCollection()
            if start not in alt_check[identifier]:
                alternatives[identifier].add_orf(orf)
                alt_check[identifier].append(start)
        print('done fetching')
        return alternatives

    def __define_identifier(self, orf_end, transcript_name=None):
        if self.subset == "Genome":
            identifier = orf_end
        else:
            identifier = f'{transcript_name}_{orf_end}'
        return identifier

    def sort_by_coordinates(self):
        """ Sorts the alternative ORFs inside self.alternatives by their start codons. """
        print('sorting by coordinates')
        for alts in self.alternatives:
            starts = []
            for alt in self.alternatives[alts]:
                if len(starts) == 0:
                    starts.append(alt)
                else:
                    if alt.strand == 'forward':
                        if starts[0].start > alt.start:
                            starts.insert(0, alt)
                        else:
                            starts.append(alt)
                    else:
                        if starts[0].start < alt.start:
                            starts.insert(0, alt)
                        else:
                            starts.append(alt)
            self.alternatives[alts] = ORFCollection()
            self.alternatives[alts].add_orfs(starts)
        # total = []
        # for stop in self.alternatives:
        #     for i in self.alternatives[stop]:
        #         total.append(f'{i.name} {i.MSPeptides} {i.start_codon} {i.proteinSequence}\n')
        # with open('teste_peptide_sort.txt', 'w') as out:
        #     out.writelines(total)

    def sort_by_atg(self):
        """ Gives priority to ATG when sorting ORFs. """
        print('sorting by atgs')
        for alts in self.alternatives:
            # print(f'Alts: {alts}')
            atgs = []
            for alt in self.alternatives[alts]:
                if len(atgs) == 0:
                    atgs.append(alt)
                else:
                    if alt.start_codon == "ATG":
                        if atgs[0].start_codon != "ATG":  # if the alt start at position 0 is not ATG, changes priority
                            atgs.insert(0, alt)  # to this ORF in the collection
                        else:  # if it is an ATG, check which one is closer to the closest upstream STOP codon
                            if alt.strand == 'forward':
                                if atgs[0].start > alt.start:
                                    atgs.insert(0, alt)
                                else:
                                    atgs.append(alt)
                            elif alt.strand == 'reverse':
                                if atgs[0].start < alt.start:
                                    atgs.insert(0, alt)
                                else:
                                    atgs.append(alt)
                    else:       # checks priority for alternative START codons. In the absence of an ATG, checks which
                        if atgs[0].start_codon != "ATG":  # alternative is closer to the closest upstream STOP codon
                            if alt.strand == 'forward':
                                if atgs[0].start > alt.start:
                                    atgs.insert(0, alt)
                                else:
                                    atgs.append(alt)
                            elif alt.strand == 'reverse':
                                if atgs[0].start < alt.start:
                                    atgs.insert(0, alt)
                                else:
                                    atgs.append(alt)
                        else:
                            atgs.append(alt)
                # print([(alt.name, alt.start, alt.start_codon) for alt in atgs])
            self.alternatives[alts] = ORFCollection()
            self.alternatives[alts].add_orfs(atgs)
            # print([(alt.name, alt.start, alt.start_codon, alt.seq) for alt in self.alternatives[alts]])
        # total = []
        # for stop in self.alternatives:
        #     for i in self.alternatives[stop]:
        #         total.append(f'{i.name} {i.MSPeptides} {i.start_codon} {i.proteinSequence}\n')
        # with open('teste_peptide_sort.txt', 'w') as out:
        #     out.writelines(total)
            # print(len(self.alternatives))

    def __fetch_codons(self, orf):
        """ :returns the nucleotide sequence of the start codon for a given ORF. """
        # print('fetch_codons')
        nucs = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        # if self.subset == "Genome":
        #     sequence = self.genome_seq[0]
        # else:
        #     sequence = orf.transcript
        sequence = self.__check_subset(orf)
        if orf.strand == 'forward':
            # if self.subset == "Genome":
            s_codon = sequence[orf.start-1: orf.start+2]
            # else:
            #     s_codon = orf.transcript[orf.start-1 : orf.start+2]
        else:
            to_reverse = sequence[orf.start-3: orf.start][::-1]
            s_codon = ""
            for nuc in to_reverse:
                s_codon += nucs[nuc]
        orf.start_codon = s_codon
        return orf

    def __check_subset(self, orf):
        if self.subset == "Genome":
            sequence = self.genome_seq[0]
        else:
            sequence = orf.transcript
        return sequence

    def __check_transcript(self, orf):
        if self.subset == "Genome":
            transcript = None
        else:
            transcript = orf.transcript
        return transcript

    def extend_orfs(self, args):
        """ note: call this function first. From now on, it's pure black magic. As this is getting kinda complex,
        to hell with python PEPs. I must remind myself to improve the readability of this chaotic mess. """
        print("extending ORFS")
        new_alts = {}
        # rev_genome = self.genome_seq[0][::-1]
        # print("extend_orfs")
        for alts in self.alternatives:
            for alt in self.alternatives[alts]:
                # print(alt.strand)
                nucs = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
                sequence = self.__check_subset(alt)
                transcript = self.__check_transcript(alt)
                if alt.strand == 'forward':
                    i = 3
                    extend = True
                    extended = 'same'
                    position = 0
                    seq = "dummy"
                    while extend and (alt.start - i) > 0:
                        s_codon = sequence[alt.start - i - 1: alt.start - i + 2]
                        # print(s_codon)
                        real_start = sequence[alt.start + 2 - i: alt.start - i + 5]  # as we stop the loop
                        # at the STOP codon, the real START codon should be 3 nucleotides downstream from that STOP
                        if real_start in args.starts:
                            extended = real_start
                            position = alt.start - i
                            # print('transcript:', alt.transcriptName, transcript)

                            seq = sequence[alt.start-i+2: alt.end]
                            # print('seq:', seq)
                            # print('forward_seq')
                            # new_alts  = self.__check_length(seq, alt, new_alts, i)
                            self.__add_extended(new_alts, position, alt, extended, seq, transcript=transcript,
                                                transcript_name=alt.transcriptName)

                        if s_codon in args.stops:
                            extend = False

                            if extended != 'same':

                                # new_alts =
                                ...
                            # print("STOP RIGHT THERE MATE")
                            # print(extended)
                        i += 3
                else:
                    i = 3
                    extend = True
                    transcript = None
                    while extend:
                        ex_start = self.genome_seq[0][alt.end - 1: alt.start + i][::-1]
                        ex_seq = self.complement(ex_start)
                        # print('rev_seq')

                        ex_codon = self.genome_seq[0][alt.start + i - 3:alt.start + i][::-1]
                        ex_codon = self.complement(ex_codon)
                        # print(ex_start)
                        i += 3
                        # print(f'codon: {ex_codon}')
                        if ex_codon in args.starts:
                            # new_alts = self.__check_length(alt=alt, new_alts=new_alts, i=i, seq=ex_seq)
                            self.__add_extended(alt=alt, new_alts=new_alts, s_codon=ex_codon, seq=ex_seq,
                                                start_pos=alt.start +i-6, transcript=transcript)
                            # print(alt.start+i-3, alt.end)
                            # print(ex_seq)
                        if ex_codon in args.stops:
                            extend = False

                    # i = 6
                    # extend = True
                    # extended = 'same'
                    # position = 0
                    # start_pos = 0
                    # end_pos = 0
                    # seq = 'dummy'
                    # print(alt.start, alt.end)
                    # while extend:
                    #     to_rev = self.genome_seq[0][alt.start - i: alt.start - i + 3][::-1]
                    #     s_codon = ""
                    #     for nuc in to_rev:
                    #         s_codon += nucs[nuc]
                    #     real_start_to_rev = self.genome_seq[0][alt.start - i + 3: alt.start + 6 - i][::-1]
                    #     real_start = ""
                    #     for nuc in real_start_to_rev:
                    #         real_start += nucs[nuc]
                    #     print(real_start)
                    #
                    #     if real_start in args.starts:
                    #         extended = real_start
                    #         position = alt.start - i
                    #         extended = real_start
                    #         position = alt.start - i
                    #         start_pos = alt.end - i + 5
                    #         end_pos = alt.end+i - 4
                    #         seq_to_rev = self.genome_seq[0][alt.end - i: alt.start][::-1]
                    #         # seq_to_rev = seq_to_rev[::-1]
                    #         seq = ""
                    #         for nuc in seq_to_rev:
                    #             seq += nucs[nuc]
                    #         print(seq)
                    #         new_alts, extend = self.__check_length(seq, alt, new_alts, i)
                    #         self.__add_extended(new_alts, position, alt, extended, seq)
                    #     if real_start in args.starts:
                    #         ...
                    #     if s_codon in args.stops:
                    #         print(f'stop: {s_codon}')
                    #         extend = False
                    #         print(f'real start: {real_start}')
                    #         if extended != 'same':
                    #             new_alts = self.__add_extended(new_alts, position, alt, extended, seq)
                    #         # print(seq)
                    #         print(alt.end + i - 4, alt.start - i + 4)
                    #         print("STOP RIGHT THERE MATE")
                    #         print(extended)
                    #     i += 3
        self.__add_new_alts(new_alts)
        print('done extending')

    @staticmethod
    def complement(seq):
        nucs = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        c = ""
        for nuc in seq:
            c += nucs[nuc]
        return c

    def __add_new_alts(self, new_alts):
        """ adds the extended ORFS to self.alternatives, so they may be compared to see which one has priority based
        on start codons composition and position. """
        for alt in new_alts:
            # print(new_alts[alt])
            # print(alt)
            # print(self.alternatives)
            # if 'Discard' not in new_alts[alt]:
            # print([(orf.name, orf.start, orf.end, orf.seq) for orf in new_alts[alt]])
            self.alternatives[str(alt)].add_orfs(new_alts[alt])
        return self

    def __check_length(self, seq, alt, new_alts, i):
        """ Deprecated """
        # extend = True
        if i > self.maxsize:
            # extend = False
            orf = ORF(name='Discard', end=alt.end, start=alt.start, seq=seq, strand=alt.strand)
            orf.shineDalgarno = alt.shineDalgarno
            orf.freeEnergy = alt.freeEnergy
            if alt.end not in new_alts:
                new_alts[alt.end] = [orf]
            else:
                new_alts[alt.end].insert(0, orf)
        return new_alts

    def __add_extended(self, new_alts, start_pos, alt, s_codon, seq, transcript, transcript_name=None):
        orf = ORF(name=f'{alt.name[:5]}_extended_{start_pos+3}-{alt.end}_{alt.strand}',
                  strand=alt.strand, start=start_pos+3, end=alt.end, seq=seq, transcript=transcript)
        orf.start_codon = s_codon
        orf.MSPeptides = alt.MSPeptides
        orf.transcriptName = transcript_name
        identifier = self.__define_identifier(orf_end=orf.end, transcript_name=transcript_name)
        if identifier not in new_alts:
            new_alts[identifier] = [orf]
        else:
            new_alts[identifier].append(orf)
        return new_alts

    def sort_by_shine(self):
        """
        Sorts the ORFs by their free energy when binding to the 16S rRNA. If no alternative START of a given STOP has
        any Shine-Dalgarno sequences (free energy >= 0) upstream from their START codon, these ORFs are not sorted at
        all.
        :return: dictionary containing ORFCollections with ORFs sorted by the presence or absence of Shine-Dalgarno
        sequences upstream from their START.
        """
        print('sorting by RBS')
        alternatives = {}
        for stop in self.alternatives:
            ignore = False
            sorted_orfs = [orf for orf in self.alternatives[stop]]
            no_sd = 0  # marks the number of sequences without a SD sequence. If this is equal to the number of
            # alternative starts, change ignore to True. Then, these ORFs should not be sorted by SD.
            for alt in self.alternatives[stop]:
                if alt.freeEnergy < 0:
                    if len(sorted_orfs) > 0:
                        if alt.freeEnergy < sorted_orfs[0].freeEnergy:
                            sorted_orfs.insert(0, alt)
                    else:
                        sorted_orfs.append(alt)
                else:
                    sorted_orfs.append(alt)
                    no_sd += 1
                if no_sd == len(self.alternatives[stop]):
                    ignore = True
            if stop not in alternatives:
                if not ignore:
                    alternatives[stop] = ORFCollection()
                    alternatives[stop].add_orfs(sorted_orfs)
                else:
                    alternatives[stop] = self.alternatives[stop]
            else:
                if not ignore:
                    alternatives[stop].add_orfs(sorted_orfs)
                else:
                    alternatives[stop] = self.alternatives[stop]
        # for stop in alternatives:
            # for alt in alternatives[stop]:
            # print('shine:', [(stop, alt.freeEnergy, alt.name, alt.start, alt.shineDalgarno) for alt in alternatives[stop]])
        # total = []
        # for stop in alternatives:
        #     i = 0
        #     for i in alternatives[stop]:
        #         if i == 0:
        #             total.append('__________________________________________\n')
        #         total.append(f'{i.name} {i.MSPeptides} {i.start_codon} {i.proteinSequence} {i.shineDalgarno} '
        #                      f'{i.freeEnergy}\n')
        #     total.append('__________________________________________\n')
        # with open('teste_peptide_sort.txt', 'w') as out:
        #     out.writelines(total)
        self.alternatives = alternatives
        return alternatives


    def __extract_peptides(self):
        peptides = self.df["Peptide"].tolist()
        # print(peptides)
        def format_pep(pep):
            fixed = ""
            # print(pep)
            for i in pep:
                if i.isalpha():
                    fixed += i
            return fixed

        fixed_peps = []
        for pep in peptides:
            fixed = format_pep(pep)
            fixed_peps.append(fixed)
        alts_with_peps = {}
        for stop in self.alternatives:
            for alt in self.alternatives[stop]:
                for pep in fixed_peps:
                    if pep not in alt.MSPeptides and pep in alt.proteinSequence:
                        # if pep in alt.proteinSequence:
                        alt.MSPeptides.append(pep)
                            # print(alt.MSPeptides, alt.proteinSequence, pep)
                            # break
                if str(stop) not in alts_with_peps:
                    alts_with_peps[str(stop)] = ORFCollection()
                    alts_with_peps[str(stop)].add_orf(alt)
                elif str(stop) in alts_with_peps:
                    alts_with_peps[str(stop)].add_orf(alt)
        return alts_with_peps

    def sort_by_peptides(self):
        print('sorting by ms peptides')
        pep_sorted = {}
        total = []
        for stop in self.alternatives:
            orfs = []
            for alt in self.alternatives[stop]:
                to_trans = Translator(alt.seq)
                alt.proteinSequence = to_trans.translate()
                alt.filter_peptides()
                if len(orfs) > 0:
                    if len(alt.MSPeptides) > len(orfs[0].MSPeptides):
                        orfs.insert(0, alt)
                    else:
                        orfs.append(alt)
                else:
                    orfs.append(alt)
            if stop not in pep_sorted:
                pep_sorted[stop] = ORFCollection()
                pep_sorted[stop].add_orfs(orfs)
            else:
                pep_sorted[stop].add_orfs(orfs)
        """ FOR DEBUGGING """
        for stop in pep_sorted:
            i = 0
            for i in pep_sorted[stop]:
                if i == 0:
                    total.append('__________________________________________\n')
                total.append(f'{i.name} {i.MSPeptides} {i.start_codon} {i.proteinSequence} {i.shineDalgarno} '
                             f'{i.freeEnergy}\n')

            total.append('__________________________________________\n')

        with open('teste_peptide_sort.txt', 'w') as out:
            out.writelines(total)
            # print([('start', orf.name, orf.MSPeptides) for orf in self.alternatives[stop]])

        self.alternatives = pep_sorted
        return pep_sorted

    def get_priorities(self):
        """ :returns the entries that have top priority for their START codons. """
        pep_list = {}
        for stop in self.alternatives:
            # pep_list.append(self.alternatives[stop])
            for orf in self.alternatives[stop]:
                pep_list[stop] = [orf]
                break
        return pep_list
