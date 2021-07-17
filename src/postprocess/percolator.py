import os
import sys
from Bio import SeqIO


class Decoy(object):
    def __init__(self, db, db_type):
        self.df = db
        self.seqs, self.entries = self.__get_seqs()
        self.reversed = []
        self.type = db_type
        self.__create_dir()
        self.path = sys.path[0]

    def __create_dir(self):
        if not os.path.exists(f"{self.type}/Percolator"):
            cmd_dir = f'mkdir {self.type}/Percolator'
            os.system(cmd_dir)

    def __get_seqs(self):
        seqs = []
        entries = []
        records = SeqIO.parse(self.df, 'fasta')
        for record in records:
            seqs.append(record.seq)
            entries.append(str(record.description))
        return seqs, entries

    def reverse_sequences(self):
        """ Reverses the amino acid sequence of a protein, except for the aa in the c-terminal. """
        for seq in self.seqs:
            # cterminus = seq[len(seq)-1]
            to_reverse = seq[:-1]
            reversed = to_reverse[::-1]
            # reversed += cterminus
            self.reversed.append(reversed)
        return self

    def add_contaminants(self):
        seqs = []
        records = SeqIO.parse(f'{self.path}/seqlib/contaminants.txt', 'fasta')
        for record in records:
            seqs.append(f">{record.id}\n{record.seq}\n")
        return seqs

    def to_fasta(self):
        out = []
        for i in range(len(self.reversed)):
            string = f">decoy_{self.entries[i]}\n{self.reversed[i]}\n"
            out.append(string)
        seqs = self.add_contaminants()
        for seq in seqs:
            to_add = seq.replace("B", "")
            to_add = to_add.replace("X", "")
            to_add = to_add.replace("Z", "")
            out.append(to_add)
        with open(f"{self.type}/Percolator/{self.type}_decoy.fasta", 'w') as fa:
            fa.writelines(out)
        return self


class PercolatorProcessing(object):
    def __init__(self, folder, filetype, args=None):
        self.folder = folder
        self.args = args
        self.target = None
        self.decoy = None
        self.filetype = filetype

    def create_metafiles(self):
        """ Creates a metafile containing the paths to the mzid files, created during the peptide search step. """
        files = os.listdir(self.folder)
        real = []
        decoy = []
        for file in files:
            if "mzid" in file:
                if "decoy" in file:
                    decoy.append(os.path.abspath(f"{self.folder}/{file}\n"))
                else:
                    real.append(os.path.abspath(f"{self.folder}/{file}\n"))
        with open(f"{self.folder}/{self.folder}_target_metafile.txt", 'w') as out:
            out.writelines(real)
        with open(f"{self.folder}/{self.folder}_decoy_metafile.txt", 'w') as decoy_out:
            decoy_out.writelines(decoy)
        self.target = f"{self.folder}/{self.folder}_target_metafile.txt"
        self.decoy = f"{self.folder}/{self.folder}_decoy_metafile.txt"
        return self

    def convert_to_pin(self):
        """ Converts the mzid files to a percolator input file using the metafiles created with create_metafiles()
        method. """
        msgf2pin_path = "msgf2pin"
        if self.args is not None:
            if self.args.msgf2pin_path is not None:
                msgf2pin_path = self.args.msgf2pin_path
        if not os.path.exists(f"{self.folder}/Percolator"):
            cmd_dir = f'mkdir {self.folder}/Percolator'
            os.system(cmd_dir)
        cmd_pin = f'{msgf2pin_path} {self.target} {self.decoy} -o {self.folder}/Percolator/{self.folder}_pin.txt' \
                  f' -F {self.filetype}_database.fasta,{self.folder}/Percolator/{self.folder}_decoy.fasta -c 2'
        os.system(cmd_pin)
        return self

    def percolate(self):
        """ Runs percolator on the converted mzid files. """
        perc_path = "percolator"
        enz = "trypsinp"
        folder = f"{self.folder}/Percolator"
        if self.args is not None:
            if self.args.enzyme is not None:
                enz = self.args.enzyme
            if self.args.percolator_path is not None:
                perc_path = self.args.percolator_path
        cmd_perc = f'{perc_path} -X {self.folder}/Percolator/{self.filetype}_percolator_out.xml --tab-out ' \
                   f'{folder}/comp_features.txt -w {folder}/feature_weights.txt -v 3 -r {folder}/pep_results.txt -m ' \
                   f'{folder}/results_psm.txt --search-input separate --results-proteins {folder}/protein_results.txt' \
                   f' --protein-enzyme {enz} --protein-report-duplicates --protein-report-fragments ' \
                   f'--spectral-counting-fdr 0.01 --picked-protein auto {folder}/{self.folder}_pin.txt'
        os.system(cmd_perc)
        return self

# os.chdir("/home/eduardo/Documents/R/R_for_proteomics/20190503_PePTID_20190530_INhAMAIS1_10uL_120min_1/")
# perc = PercolatorProcessing("Genome", filetype="genome")
# perc.create_metafiles().convert_to_pin()
# perc.percolate()
# genome_decoy = Decoy(db="genome_ORFs.fasta", db_type="Genome")
# genome_decoy.reverse_sequences().to_fasta()
