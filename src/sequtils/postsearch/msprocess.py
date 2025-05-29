# Copyright © 2021-2025 Eduardo Vieira de Souza
# Copyright © 2021-2025 Adriana Canedo
# Copyright © 2021-2025 Cristiano Valim Bizarro
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

from ..__helpers import EngineError, PercolatorProteinsError, PercolatorPSMError, UProteinsError
from ..orflib import ORF, ORFCollection
from .peplib import Peptide, PeptideCollection
from ..utilities import find_coords, GFFReader, findnth


class PercolatorUP(object):
    def __init__(self, proteins=None, peptides=None):
        self.proteins = proteins
        self.peptides = peptides
        self.__check_files()
        self.proteinsDataFrame = pd.read_csv(self.proteins, sep="\t")
        self.PSMDataFrame = pd.read_csv(self.peptides, sep="\t", usecols=['PSMId', 'score',	'q-value', 'posterior_error_prob', 'peptide'])
        self.__fix_columns()

    def __check_files(self):
        if self.peptides is None:
            raise PercolatorPSMError

    def __fix_columns(self):
        """ Puts all protein Ids of percolator output in a single column. """
        ids = []
        with open(self.peptides, 'r') as psm:
            lines = psm.readlines()
            for i in range(len(lines)):
                if i > 0:
                    proteins = lines[i].split("\t")[5:]
                    sep = ","
                    proteins = sep.join(proteins)
                    proteins = proteins.rstrip()
                    ids.append(proteins)
        self.PSMDataFrame.insert(5, "proteinIds", ids)

    def filter_proteins(self, **kwargs):
        """ Filter the protein results by the specified cutoffs. Accepted arguments are 'q-value', 'eprob', 'spec',
        and 'spec-all'. """
        if self.proteins is None:
            raise PercolatorProteinsError
        qvalue = kwargs.get("q-value")
        eprob = kwargs.get("eprob")
        spec = kwargs.get("spec")
        spec_all = kwargs.get("spec-all")
        df = self.proteinsDataFrame
        if qvalue is not None:
            df = df[df["q-value"] <= qvalue]
        if eprob is not None:
            df = df[df["posterior_error_prob"] <= eprob]
        if spec is not None:
            df = df[df["spec_count_unique"] >= spec]
        if spec_all is not None:
            df = df[df["spec_all"] >= spec_all]
        self.proteinsDataFrame = df
        return self

    def filter_psm(self, **kwargs):
        """ Filters the PSM results by the specified cutoffs. Accepted arguments are 'score', 'q-value' and 'eprob. """
        score = kwargs.get("score")
        qvalue = kwargs.get("q-value")
        eprob = kwargs.get("eprob")
        df = self.PSMDataFrame
        if score is not None:
            df = df[df['score'] <= score]
        if qvalue is not None:
            df = self.peptides[self.peptides['q-value'] <= qvalue]
        if eprob is not None:
            df = self.peptides[self.peptides['posterior_error_prob'] <= eprob]
        self.PSMDataFrame = df
        return self

    def filter_transcriptome_unique(self, method="strict", fasta=None, keep='unique', stringtie=None, orf_db=None):
        """ Filter the unique peptides in the percolator output. Different methods are available: 'strict' removes any
        PSM that matches two or more proteins. 'uproteins' uses uProteInS method of classifying unique peptides.
        'alternative' also classifies as unique any PSM that matches more than two proteins at the same genomic loci,
        as long as they have distinct start codons. Both 'uproteins' and 'alternative' methods requires a genbank file
        or a tab-separated data frame containing coordinates of all in the set. 'keep' accepts either 'unique' or
        'all'. If set to unique, it removes all non-unique peptides. Otherwise, it still checks whether the peptide
        is unique or not, but keeps this information in the data frame regardless."""
        ids = self.PSMDataFrame["proteinIds"].tolist()
        peptides = self.PSMDataFrame["peptide"].tolist()
        pep_dict = {}
        if method == 'strict':
            checked_peps = self.__strict(ids, peptides)
        elif method == 'uproteins':
            checked_peps = self.__uproteins(ids, peptides, stringtie, orf_db)
        for pep in checked_peps:
            pep_dict[pep.seq] = pep
        unique = []
        for i in range(len(peptides)):
            if pep_dict[peptides[i]].unique:
                unique.append(True)
            else:
                unique.append(False)
        self.PSMDataFrame.insert(6, "Unique Peptide", unique)
        if keep == 'unique':
            self.PSMDataFrame = self.PSMDataFrame[self.PSMDataFrame["Unique Peptide"] == True]
        return self

    def filter_genome_unique(self, method='strict', gff=None, orf_db=None, keep='all'):
        ids = self.PSMDataFrame["proteinIds"].tolist()
        peptides = self.PSMDataFrame["peptide"].tolist()
        pep_dict = {}
        if method == 'uproteins':
            checked_peps = self.__uproteins_dna(ids, peptides, orf_db, gff)
        for pep in checked_peps:
            pep_dict[pep.seq] = pep
        unique = []
        for i in range(len(peptides)):
            if pep_dict[peptides[i]].unique:
                unique.append(True)
            else:
                unique.append(False)
        self.PSMDataFrame.insert(6, "Unique Peptide", unique)
        if keep == 'unique':
            self.PSMDataFrame = self.PSMDataFrame[self.PSMDataFrame["Unique Peptide"] == True]
        return self

    def __uproteins_dna(self, proteins, peptides, orf_db, gff):
        genes = GenomeMiner(orf_db=orf_db, gff=gff)
        novel_orfs = genes.get_orfs(subset='novel')
        ref_orfs = genes.get_orfs(subset='refseq')
        novel_dict = self.__get_orf_dict(novel_orfs)
        ref_dict = self.__get_orf_dict(ref_orfs)
        pep_list = []
        for pro, pep in zip(proteins, peptides):
            peptide = Peptide(sequence=pep)
            for protein in pro.split(","):
                if protein in novel_dict:
                    start = novel_dict[protein].start
                    end = novel_dict[protein].end
                    peptide.add_novel_spec(start, end)
                elif protein in ref_dict:
                    start = ref_dict[protein].start
                    end = ref_dict[protein].end
                    peptide.add_ref_spec(start, end)
            # peptide.add_novel_spec()
            peptide.check_loci()
            pep_list.append(peptide)
        return PeptideCollection(pep_list)


    def get_ups(self, method='strict', stringtie=None, orf_db=None):
        """ Filters unique peptides just as the filter_unique() method does, but returns an instance of ORFCollection
        containing information about the peptides of each ORF instead of saving this information in a data frame."""
        ids = self.PSMDataFrame["proteinIds"].tolist()
        peptides = self.PSMDataFrame["peptide"].tolist()
        if method == 'strict':
            checked_peps = self.__strict(ids, peptides)
        elif method == 'uproteins':
            checked_peps = self.__uproteins(ids, peptides, stringtie, orf_db)
        return checked_peps

    def __strict(self, proteins, peptides):
        pep_list = []
        for pep, protein in zip(peptides, proteins):
            peptide = Peptide(pep)
            if len(protein.split(",")) > 1:
                peptide.unique = False
            else:
                peptide.unique = True
            pep_list.append(peptide)
        return PeptideCollection(pep_list)

    def __uproteins(self, proteins, peptides, stringtie, orf_db):
        genes = TranscriptomeMiner(orf_db=orf_db, gff=stringtie)
        novel_orfs = genes.get_orfs(subset='novel')
        ref_orfs = genes.get_orfs(subset='refseq')
        novel_dict = self.__get_orf_dict(novel_orfs)
        ref_dict = self.__get_orf_dict(ref_orfs)
        pep_list = []
        for pro, pep in zip(proteins, peptides):
            peptide = Peptide(sequence=pep)
            for protein in pro.split(","):
                if protein in novel_dict:
                    start = novel_dict[protein].start
                    end = novel_dict[protein].end
                    peptide.add_novel_spec(start, end)
                elif protein in ref_dict:
                    start = ref_dict[protein].start
                    end = ref_dict[protein].end
                    peptide.add_ref_spec(start, end)
            # peptide.add_novel_spec()
            peptide.check_loci()
            pep_list.append(peptide)
        return PeptideCollection(pep_list)

    def __get_orf_dict(self, orfset):
        orf_dict = {}
        for orf in orfset:
            orf_dict[orf.name] = orf
        return orf_dict

    def save_psms(self, filename="psms.txt"):
        self.PSMDataFrame.to_csv(filename, sep="\t", index=False)


class RNAMiner(object):
    """ Extracts information about the predicted ORFs locus in the genome using the transcript information."""
    def __init__(self, orf_db=None, genbank=None, transcripts=None):
        """ transcripts accepts the output GTF file from stringtie. If it is not specified, this method will assume
        the ORFs were predicted from the genome. """
        self.db = orf_db
        self.gb = genbank
        self.transcripts = transcripts
        # self.__check_files()

    def __check_files(self):
        if self.db is None or self.gb is None:
            raise UProteinsError

    def _get_orfs(self):
        orf_list = []
        records = SeqIO.parse(self.db, 'fasta')
        for record in records:
            if 'uproteins' in record.id:

                start, end = find_coords(str(record.description))
                orf = ORF(name=record.description, seq=record.seq, start=start, end=end)
                # no_orf = record.id[4:]
                # pos = no_orf.find("_")
                # trans = no_orf[:pos]
                # end = findnth(trans, ".", 2)
                # orf.transcript = trans[:end-1]
                # orf_dict[orf.name] = orf
                gene = record.id
                if 'uproteins' not in gene:
                    pos = gene[:-1].rfind("_")
                    new = gene[:pos]
                    pos2 = new.rfind("_")
                    new = new[:pos2]
                    orf_name = new[9:]
                else:
                    no_orf = gene[4:]
                    pos = no_orf.find("_")
                    trans = no_orf[:pos]
                    end = findnth(trans, ".", 2)
                    orf_name = trans[:end - 1]
                orf.transcript = orf_name
                orf_list.append(orf)
        return orf_list

    def genome_locus(self, subset):
        """ Subset is either 'novel' or 'refseq'. Pattern refers to the pattern present in the fasta file that
        informs if a ORF is predicted or not. """
        rna = GFFReader(gff=self.transcripts)
        # pd.set_option('max_rows', 50)
        novel, ref = rna.find_novel()
        db = self._get_orfs()
        located_orfs = []
        subset_dict = {'novel': novel, 'refseq': ref}
        df = subset_dict[subset]
        records = SeqIO.parse(self.db, 'fasta')
        if subset == 'novel':
            for record in records:
                if 'uproteins' in record.id:
                    df = df[df["name"] == db[orf].transcript]
                    # position of the transcript in the genome
                    rna_genome_start = df["start"].tolist()[0]

                    rna_genome_end = df["end"].tolist()[0]

                    # position of the ORF in the genome
                    orf_genome_start = int(rna_genome_start) + int(db[orf].start)-1
                    orf_genome_end = int(rna_genome_start) + int(db[orf].end)-1

                    # strand info
                    rna_strand = df["strand"].tolist()[0]
                    if rna_strand == '-':
                        strand = 'reverse'
                    else:
                        strand = 'forward'

                    # adds info to a new ORF object
                    corrected = ORF(name=orf, transcript=db[orf].transcript, start=orf_genome_start, end=orf_genome_end,
                                    strand=strand, seq=db[orf].seq)
                    located_orfs.append(corrected)
            return ORFCollection().add_orfs(located_orfs)
        else:
            for orf in db:
                if '|' in orf:
                    df = df[df["name"] == db[orf].transcript]
                    rna_genome_start = df["start"].tolist()[0]
                    rna_genome_end = df["end"].tolist()[0]

                    # position of the ORF in the genome
                    orf_genome_start = int(rna_genome_start) + int(db[orf].start) - 1
                    orf_genome_end = int(rna_genome_start) + int(db[orf].end) - 1

                    # strand info
                    rna_strand = df["strand"].tolist()[0]
                    if rna_strand == '-':
                        strand = 'reverse'
                    else:
                        strand = 'forward'

                    # adds info to a new ORF object
                    corrected = ORF(name=orf, transcript=db[orf].transcript, start=orf_genome_start, end=orf_genome_end,
                                    strand=strand, seq=db[orf].seq)
                    located_orfs.append(corrected)
            return ORFCollection().add_orfs(located_orfs)


    def __check_(self, orf, start, end):
        if self.transcripts is None:
            if 'reverse' in orf:
                strand = 'reverse'
            else:
                strand = 'forward'
            return strand
        else:
            pass


class TranscriptomeMiner(object):
    def __init__(self, orf_db=None, gff=None):
        """ 'gff' accepts a stringtie GFF3 file. """
        self.db = orf_db
        self.gff = gff
        self._parse_ref()

    def _parse_ref(self):
        anno = GFFReader(refseq=self.gff)
        orfs = anno.find_annotated_rna()
        self.annotated = orfs

    def find_gene_name(self, gene):
        if 'uproteins' not in gene:
            pos = gene[:-1].rfind("_")
            new = gene[:pos]
            pos2 = new.rfind("_")
            new = new[:pos2]
            return new[9:]
        else:
            no_orf = gene[4:]
            pos = no_orf.find("_")
            trans = no_orf[:pos]
            end = findnth(trans, ".", 2)
            return no_orf[:end-1]

    def get_orfs(self, subset='novel'):
        records = SeqIO.parse(self.db, 'fasta')
        orf_list = []
        if subset == 'novel':
            for record in records:
                    if 'ORF' in record.id:
                        end_coords = record.id.rfind("_")
                        id_to_end = record.id[:end_coords]
                        start_coords = id_to_end.rfind("_") + 1
                        coords = id_to_end[start_coords:].split("-")

                        if 'forward' in record.id:
                            start = coords[0]
                            end = coords[1]
                            strand = 'forward'
                        else:
                            start = coords[1]
                            end = coords[0]
                            strand = 'reverse'

                        orf = ORF(name=record.id, seq=record.seq, start=start, end=end, strand=strand)
                        orf_list.append(orf)
            return ORFCollection().add_orfs(orf_list)
        else:
            for record in records:
                if "MSMEG" in record.id:
                    df = self.annotated[self.annotated["name"] == record.id]
                    start = df["start"].tolist()[0]
                    end = df["end"].tolist()[0]
                    strand = df["strand"].tolist()[0]
                    if strand == '-':
                        strand = 'reverse'
                    elif strand == '+':
                        strand = 'forward'
                    orf = ORF(name=record.id, seq=record.seq, start=start, end=end, strand=strand)
                    orf_list.append(orf)
            return ORFCollection().add_orfs(orf_list)



class GenomeMiner(object):
    def __init__(self, orf_db=None, gff=None):
        """ 'gff' accepts a RefSeq GFF3 file. """
        self.db = orf_db
        self.gff = gff
        self._parse_ref()

    def _parse_ref(self):
        anno = GFFReader(refseq=self.gff)
        orfs = anno.find_annotated()
        self.annotated = orfs

    def find_gene_name(self, gene):
        if 'uproteins' not in gene:
            pos = gene[:-1].rfind("_")
            new = gene[:pos]
            pos2 = new.rfind("_")
            new = new[:pos2]
            return new[9:]
        else:
            no_orf = gene[4:]
            pos = no_orf.find("_")
            trans = no_orf[:pos]
            end = findnth(trans, ".", 2)
            return no_orf[:end-1]

    def get_orfs(self, subset='novel'):
        records = SeqIO.parse(self.db, 'fasta')
        orf_list = []
        if subset == 'novel':
            for record in records:
                    if 'ORF' in record.id:
                        end_coords = record.id.rfind("_")
                        id_to_end = record.id[:end_coords]
                        start_coords = id_to_end.rfind("_") + 1
                        coords = id_to_end[start_coords:].split("-")

                        if 'forward' in record.id:
                            start = coords[0]
                            end = coords[1]
                            strand = 'forward'
                        else:
                            start = coords[1]
                            end = coords[0]
                            strand = 'reverse'

                        orf = ORF(name=record.id, seq=record.seq, start=start, end=end, strand=strand)
                        orf_list.append(orf)
            return ORFCollection().add_orfs(orf_list)
        else:
            for record in records:
                if "MSMEG" in record.id:
                    df = self.annotated[self.annotated["name"] == record.id]
                    start = df["start"].tolist()[0]
                    end = df["end"].tolist()[0]
                    strand = df["strand"].tolist()[0]
                    if strand == '-':
                        strand = 'reverse'
                    elif strand == '+':
                        strand = 'forward'
                    orf = ORF(name=record.id, seq=record.seq, start=start, end=end, strand=strand)
                    orf_list.append(orf)
            return ORFCollection().add_orfs(orf_list)


class UProteins(object):
    def __init__(self):
        pass