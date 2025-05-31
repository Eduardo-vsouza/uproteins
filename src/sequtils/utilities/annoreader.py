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


class GFFReader(object):
    def __init__(self, orf_db=None, refseq=None, gff=None):
        """ transcripts accepts the output GTF file from stringtie. If it is not specified, this method will assume
        the ORFs were predicted from the genome. 'refseq' accepts a RefSeq GFF file. """
        self.db = orf_db
        self.refseq = refseq
        self.gff = gff

    def find_novel(self):
        """ Look up transcript information in the StringTie GTF output file. Then, it separates them into novel and
         already-annotated transcripts. Returns two pandas data frames. """
        if self.gff is not None:
            # pd.set_option('max_colwidth', 300)
            df = pd.read_csv(self.gff, sep='\t', header=2)

            df.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
            novel = df[df['attributes'].str.contains('ref_gene_id ') == False]
            novel = novel[novel['source'] == 'StringTie']
            ref = df[(df['attributes'].str.contains('ref_gene_id')) | (df['source'] == 'RefSeq')]
            novel_ids, ref_ids = self.__remove_unused_cols(novel, ref)
            novel.insert(0, "name", novel_ids)
            ref.insert(0, "name", ref_ids)
            return novel, ref
        else:
            return None

    def find_annotated(self):
        """ Look up transcript information in the RefSeq GTF file. Then, it gets the already-annotated transcripts.
         Returns a pandas data frame. """
        if self.refseq is not None:
            pd.set_option('max_columns', None)
            pd.set_option('max_colwidth', 300)
            df = pd.read_csv(self.refseq, sep='\t', header=2)

            df.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
            ids = self.__rename_anno(df)
            df.insert(0, "name", ids)
            df = df[(df["source"] == "Protein Homology") & (df["feature"] == "CDS")]
            return df

    def find_annotated_rna(self):
        """ Look up transcript information in the RefSeq GTF file. Then, it gets the already-annotated transcripts.
         Returns a pandas data frame. """
        if self.refseq is not None:
            pd.set_option('max_columns', None)
            pd.set_option('max_colwidth', 300)
            df = pd.read_csv(self.refseq, sep='\t', header=2)

            df.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
            ids = self.__rename_anno(df)
            df.insert(0, "name", ids)
            df = df[(df["source"] == "Protein Homology") & (df["feature"] == "CDS")]
            return df

    def __rename_anno(self, df):
        """ Renames annotated GFF file attributes column. """
        # df = df[(df["source"] == "RefSeq") & (df["feature"] == "gene")]
        df = df["attributes"].tolist()
        ids = []

        for i in range(len(df)):
            print(df[i])
            gene = df[i].split(";")[1][12:]
            print(gene)
            ids.append(gene)
        return ids

    def __remove_unused_cols(self, novel, ref):
        """ For StringTie GFF files. """
        cols = ['start', 'end', 'attributes']
        novel = novel[cols]
        novel_at = novel['attributes'].tolist()
        novel_ids = self.__rename_attributes(novel_at)
        ref = ref[cols]
        ref_at = ref['attributes'].tolist()
        ref_ids = self.__rename_attributes(ref_at)
        return novel_ids, ref_ids

    def __rename_attributes(self, attributes):
        ids = []
        for i in range(len(attributes)):
            name = attributes[i]
            if 'ref_gene_id' in name:
                start = name.find('ref_gene_id "')+18
            else:
                start = name.find('gene_id "')+9
            end = name[start:].find(";")-1
            gene_name = name[start:][:end]
            print(gene_name)
            ids.append(gene_name)
        return ids