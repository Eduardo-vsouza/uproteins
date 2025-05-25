import pandas as pd
from Bio import SeqIO


class FastaConverter(object):
    def __init__(self, fasta=None, conversion_file=None, gff=None):
        self.fasta = fasta
        self.cFile = conversion_file
        self.refSeq = gff
        self.gff = self.__gff_info()

    def __gff_info(self):
        """ Must specify a gff file. """
        df = pd.read_csv(self.refSeq, sep='\t', header=2)
        df.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
        df = df[df["source"] == 'Protein Homology']
        return df

    def convert_entries(self, sep='\t', **kwargs):
        """ Convert all protein entries in the percolator output to the Uniprot. If 'pattern' is specified as an
        :argument, this will check only for target IDs with this pattern during the conversion. This is useful when
        the uniprot conversor returns more than a single target for a ID. """
        if kwargs.get("pattern"):
            pattern = kwargs.get("pattern")
        else:
            pattern = ""
        id_dict = {}
        df = pd.read_csv(self.cFile, sep=sep)
        df = df[df["To"].str.contains(pattern)]
        ids = df["From"].tolist()
        target = df["To"].tolist()
        for i in range(len(ids)):
            gff = self.gff[self.gff["attributes"].str.contains(target[i])]
            if len(gff["attributes"].tolist()) > 0:
                attrs = gff["attributes"].tolist()[0]
                gene_name = attrs.split(";")
                gene = gene_name[1][12:]
                id_dict[ids[i]] = gene
        fasta_file = []
        records = SeqIO.parse(self.fasta, 'fasta')
        for record in records:
            if "|" in record.id:
                name = record.id.split("|")[2]
                if name in id_dict:
                    converted = id_dict[name]
                else:
                    converted = name
            else:
                converted = record.id
            fasta_file.append(f'>{converted}\n{record.seq}\n')
        with open('converted_genome_database.fasta', 'w') as fa:
            fa.writelines(fasta_file)



class PercolatorConverter(object):
    def __init__(self, pout=None, handle='psm', conversion_file=None, gff=None):
        """ 'pout' is a percolator output PSM or Protein table. Handle specifies whether this table is 'psm' or
        'protein'"""
        self.handle = handle
        self.pout = pout
        self.__check_handle()
        self.cFile = conversion_file
        self.refSeq = gff
        self.gff = self.__gff_info()

    def __check_handle(self):
        if self.handle == 'psm':
            self.dataFrame = pd.read_csv(self.pout, sep="\t", usecols=[0, 1, 2, 3, 4, 5])
            self.dataFrame = self.dataFrame.drop('proteinIds', axis=1)
            self.__fix_columns()
        elif self.handle =='protein':
            self.dataFrame = pd.read_csv(self.pout, sep="\t")

    def __fix_columns(self):
        """ Puts all protein Ids of percolator output in a single column. """
        ids = []
        with open(self.pout, 'r') as psm:
            lines = psm.readlines()
            for i in range(len(lines)):
                if i > 0:
                    proteins = lines[i].split("\t")[5:]
                    sep = ","
                    proteins = sep.join(proteins)
                    proteins = proteins.rstrip()
                    ids.append(proteins)
        # print(self.dataFrame)
        self.dataFrame.insert(5, "proteinIds", ids)

    def convert_entries(self, sep='\t', **kwargs):
        """ Convert all protein entries in the percolator output to the Uniprot. If 'pattern' is specified as an
        :argument, this will check only for target IDs with this pattern during the conversion. This is useful when
        the uniprot conversor returns more than a single target for an ID. """
        # if kwargs.get("pattern"):
        #     pattern = kwargs.get("pattern")
        # else:
        #     pattern = ""
        # id_dict = {}
        # df = pd.read_csv(self.cFile, sep=sep)
        # df = df[df["To"].str.contains(pattern)]
        # ids = df["From"].tolist()
        # target = df["To"].tolist()
        # for i in range(len(ids)):
        #     gff = self.gff[self.gff["attributes"].str.contains(target[i])]
        #     if len(gff["attributes"].tolist()) > 0:
        #         attrs = gff["attributes"].tolist()[0]
        #         gene_name = attrs.split(";")
        #         gene = gene_name[1][12:]
        #         id_dict[ids[i]] = gene
        # converted = self.__convert_proteins(id_dict)
        # self.dataFrame = self.dataFrame.drop('proteinIds', axis=1)
        # self.dataFrame.insert(5, "proteinIds", converted)
        self.dataFrame.to_csv(f'{kwargs.get("output")}.txt', sep="\t", index=False)

    def __convert_proteins(self, id_dict):
        proteins = self.dataFrame["proteinIds"].tolist()
        converted = []
        for i in range(len(proteins)):
            p_set = proteins[i].split(",")
            c_set = ""
            for pro in range(len(p_set)):
                if 'ORF' not in p_set[pro]:
                    if p_set[pro] != "":
                        pos = p_set[pro].rfind("|") + 1
                        name = p_set[pro][pos:]
                        if name in id_dict:
                            conv = id_dict[name]
                        else:
                            conv = p_set[pro]
                        if pro > 0:
                            c_set += f',{conv}'
                        else:
                            c_set += f'{conv}'
                else:
                    if p_set[pro] != "":
                        if pro > 0:
                            c_set += f',{p_set[pro]}'
                        else:
                            c_set += f'{p_set[pro]}'
            converted.append(c_set)
        return converted

    def __gff_info(self):
        """ Must specify a gff file. """
        df = pd.read_csv(self.refSeq, sep='\t', header=2)
        df.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
        df = df[df["source"] == 'Protein Homology']
        return df