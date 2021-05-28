import pandas as pd
from Bio import SeqIO

from . import ORF, ORFCollection
from ..conversion import Translator


class AltCodons(object):
    def __init__(self, file, genome, maxsize):
        """ I hate this code """
        self.maxSize = maxsize
        self.df = pd.read_csv(file, sep='\t')
        self.coordinates = self.df["Genome Coordinates"].tolist()
        self.names = self.df["Protein"].tolist()
        self.proteinSequences = self.df["db entry"].tolist()
        print(self.proteinSequences)
        self.__genome_records = SeqIO.parse(genome, 'fasta')
        self.genome_seq = [str(record.seq) for record in self.__genome_records]

        self.alternatives = self.__fetch_orfs()
        self.alternatives = self.__extract_peptides()
        for stop in self.alternatives:
            for alt in self.alternatives[stop]:
                print(alt.MSPeptides, alt.name)
                break

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

    def __fetch_orfs(self):
        """
        :returns a dictionary containing all ORFs with alternative START codons for a given STOP codon.
        """
        alt_check = {}
        alternatives = {}
        for i in range(len(self.names)):
            start, end, strand = self.__split_coords(i)
            seq = self.genome_seq[0][start - 1: end+1]
            orf = ORF(name=self.names[i], start=int(start), end=int(end), seq=seq, strand=strand, protein_sequence=self.proteinSequences[i])
            orf = self.__fetch_codons(orf)
            if end not in alt_check:
                alt_check[end] = []
            # if start not in alt_check[end]:
            #     alt_check[end].append(start)
            if end not in alternatives:
                alternatives[end] = ORFCollection()
            if start not in alt_check[end]:
                alternatives[end].add_orf(orf)
                alt_check[end].append(start)
        return alternatives

    def sort_by_coordinates(self):
        """ Sorts the alternative ORFs inside self.alternatives by their start codons. """
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
        nucs = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        if orf.strand == 'forward':
            s_codon = self.genome_seq[0][orf.start-1: orf.start+2]
        else:
            to_reverse = self.genome_seq[0][orf.start-3: orf.start][::-1]
            s_codon = ""
            for nuc in to_reverse:
                s_codon += nucs[nuc]
        orf.start_codon = s_codon
        return orf

    def extend_orfs(self, args):
        """ note: call this function first. From now on, it's pure black magic. As this is getting kinda complex,
        to hell with python PEPs. I must remind myself to improve the readability of this chaotic mess. """
        new_alts = {}
        rev_genome = self.genome_seq[0][::-1]
        for alts in self.alternatives:
            for alt in self.alternatives[alts]:
                # print(alt.strand)
                nucs = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
                if alt.strand == 'forward':
                    i = 6
                    extend = True
                    extended = 'same'
                    position = 0
                    seq = "dummy"
                    while extend and (alt.start - i - 1) > 0:
                        s_codon = self.genome_seq[0][alt.start - i - 1: alt.start - i + 2]
                        # print(s_codon)
                        real_start = self.genome_seq[0][alt.start + 2 - i: alt.start - i + 5]  # as we stop the loop
                        # at the STOP codon, the real START codon should be 3 nucleotides downstream from that STOP
                        if real_start in args.starts.split(","):
                            extended = real_start
                            position = alt.start - i

                            seq = self.genome_seq[0][position+2: alt.end+1]
                            # new_alts  = self.__check_length(seq, alt, new_alts, i)
                            self.__add_extended(new_alts, position, alt, extended, seq)
                        if s_codon in args.stops.split(","):
                            extend = False

                            if extended != 'same':

                                # new_alts =
                                ...
                            # print("STOP RIGHT THERE MATE")
                            # print(extended)
                        i += 3
                else:
                    i = 6
                    extend = True
                    while extend:
                        ex_start = self.genome_seq[0][alt.end - 1: alt.start + i][::-1]
                        ex_seq = self.complement(ex_start)
                        ex_codon = self.genome_seq[0][alt.start + i - 3:alt.start + i][::-1]
                        ex_codon = self.complement(ex_codon)
                        # print(ex_start)
                        i += 3
                        # print(f'codon: {ex_codon}')
                        if ex_codon in args.starts.split(","):
                            # new_alts = self.__check_length(alt=alt, new_alts=new_alts, i=i, seq=ex_seq)
                            self.__add_extended(alt=alt, new_alts=new_alts, s_codon=ex_codon, seq=ex_seq, start_pos=alt.start +i-6)
                            # print(alt.start+i-3, alt.end)
                            # print(ex_seq)
                        if ex_codon in args.stops.split(","):
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
                    #     if real_start in args.starts.split(","):
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
                    #     if real_start in args.starts.split(","):
                    #         ...
                    #     if s_codon in args.stops.split(","):
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
        # extend = True
        if i > self.maxSize:
            # extend = False
            orf = ORF(name='Discard', end=alt.end, start=alt.start, seq=seq, strand=alt.strand)
            orf.shineDalgarno = alt.shineDalgarno
            orf.freeEnergy = alt.freeEnergy
            if alt.end not in new_alts:
                new_alts[alt.end] = [orf]
            else:
                new_alts[alt.end].insert(0, orf)
        return new_alts

    @staticmethod
    def __add_extended(new_alts, start_pos, alt, s_codon, seq):
        orf = ORF(name=f'{alt.name[:5]}_extended_{start_pos+3}-{alt.end}_{alt.strand}',
                  strand=alt.strand, start=start_pos+3, end=alt.end, seq=seq)
        orf.start_codon = s_codon
        orf.MSPeptides = alt.MSPeptides
        if alt.end not in new_alts:
            new_alts[alt.end] = [orf]
        else:
            new_alts[alt.end].append(orf)
        return new_alts

    def sort_by_shine(self):
        """
        Sorts the ORFs by their free energy when binding to the 16S rRNA. If no alternative START of a given STOP has
        any Shine-Dalgarno sequences (free energy >= 0) upstream from their START codon, these ORFs are not sorted at
        all.
        :return: dictionary containing ORFCollections with ORFs sorted by the presence or absence of Shine-Dalgarno
        sequences upstream from their START.
        """
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
        pep_sorted = {}
        total = []
        for stop in self.alternatives:
            orfs = []
            for alt in self.alternatives[stop]:
                # if 'extend' not in alt.name:
                # if alt.name == 'gORF__extended_3593023-3592712_reverse':  # MUST REMOVE AFTERWARDS
                #     alt.MSPeptides = ['FAKEPEPTIDE', 'blibliblo']  # MUST DELETE
                to_trans = Translator(alt.seq)
                alt.proteinSequence = to_trans.translate()
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
        """ FOR DEBUGGING 
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
        """
        self.alternatives = pep_sorted
        return pep_sorted

    def get_priorities(self):
        """ :returns the entries that have top priority for their START codons. """
        pep_list = []
        for stop in self.alternatives:
            # pep_list.append(self.alternatives[stop])
            for orf in self.alternatives[stop]:
                pep_list.append(orf.name)
                break
        return pep_list
