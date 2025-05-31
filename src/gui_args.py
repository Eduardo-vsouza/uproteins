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


import os

from src.database import database_generator as dg
import peptide_search as ps
import results_new_approach as pms


class Arguments(object):
    def __init__(self, args_dic):
        self.dic = args_dic
        self.Transcriptome = None
        self.outdir = "/ProteInS"

    def set_assembly_standard(self):
        """Sets standard arguments for assembly mode. """
        self.rna_right = None
        self.rna_left = None
        self.rna_single = None
        self.lib_type = None
        self.k_mer = None
        self.cpu = None
        self.memory = None

    def set_db_standard(self):
        """ Sets standard arguments for database mode. """
        self.genome = ""
        self.proteome = ""
        self.starts = ""
        self.ends = ""
        self.maxsize = None
        self.maxsize = None

    def set_pepsearch_standard(self):
        self.Mass_spec = None
        self.modification1 = None
        self.modification2 = None
        self.modification3 = None
        self.modification4 = None
        self.memory = None
        self.mods = None
        self.enzyme = None
        self.tolerance = None
        self.fragmentation = None
        self.protocol = None
        self.instrument = None
        self.decoy = None
        self.termini = None
        self.isotope = None
        self.length = None
        self.charge = None
        self.matches = None
        self.numMods = 3
        self.missCleavages = None
        self.irrCleavages = None
        self.ParentError = 30
        self.Fdr = 0.01

    def set_postms_standard(self):
        self.Mass_spec = None
        self.genome = None
        self.proteome = None
        self.genbank = None


    def get_args_fullpath(self, mode):
        fullpaths = {'assembly': 'rna_left,rna_right,rna_single,outdir', 'database': 'genome,proteome,outdir',
                     'peptide_search': 'Mass_spec,outdir', 'postms': 'Mass_spec,outdir,proteome,genome,genbank'}
        for arg in self.dic:
            if arg in fullpaths[mode].split(","):
                self.dic[arg] = os.path.abspath(self.dic[arg])


    def assembly(self):
        """ Set input arguments for assembly mode. Calls run_assembly() function to perform the Transcriptome
         assembly with arguments set in this assembly() function."""
        self.set_assembly_standard()
        if "rna_left" not in self.dic and "rna_single" not in self.dic:
            print("Please inform both paired-end reads or a single-end read. Not enough reads were informed.")
        else:
            if "rna_left" in self.dic:
                self.rna_left = self.dic["rna_left"]
            if "rna_right" in self.dic:
                self.rna_right = self.dic["rna_right"]
            if "rna_single" in self.dic:
                self.rna_single = self.dic["rna_single"]
            if "lib_type" in self.dic:
                self.lib_type = self.dic["lib_type"]
            if "k_mer" in self.dic:
                self.k_mer = self.dic["k_mer"]
            if "cpu" in self.dic:
                self.cpu = self.dic["cpu"]
            if "memory" in self.dic:
                self.memory = self.dic["memory"]
            if "outdir" in self.dic:
                self.outdir = os.path.abspath(self.dic["outdir"])
            self.get_args_fullpath('assembly')
            self.run_assembly()

    def database(self):
        """ Set input argument for database mode. Calls run_db() function to generate the databases with the arguments
        set in this database() function. """
        self.set_db_standard()
        if "genome" not in self.dic or "proteome" not in self.dic:
            print("Please select both a Genome and a Proteome fasta file.")
        else:
            self.genome = self.dic["genome"]
            self.proteome = self.dic["proteome"]
            if "outdir" in self.dic:
                self.outdir = os.path.abspath("./")
            if "starts" in self.dic:
                self.starts = self.dic["starts"]
            if "ends" in self.dic:
                self.ends = self.dic["ends"]
            if "maxsize" in self.dic:
                self.maxsize = self.dic["maxsize"]
            if "minsize" in self.dic:
                self.minsize = self.dic["minsize"]
            if "Transcriptome" in self.dic:
                self.Transcriptome = self.dic["Transcriptome"]
            self.get_args_fullpath('database')
            print(self.dic)
            self.run_db()

    def peptide_search(self):
        """ Set input arguments for peptide search mode. Calls run_pep_search() function to match spectra against the
        predicted ORF sequences. """
        self.set_pepsearch_standard()
        if "Mass_spec" not in self.dic:
            print("Please inform the folder containing your raw MS data files.")
        else:
            self.Mass_spec = os.path.abspath(self.dic["Mass_spec"])
            if "outdir" in self.dic:
                self.outdir = self.dic["outdir"]
            if "modification1" in self.dic:
                self.modification1 = self.dic["modification1"]
            if "modification2" in self.dic:
                self.modification2 = self.dic["modification2"]
            if "modification3" in self.dic:
                self.modification3 = self.dic["modification3"]
            if "modification4" in self.dic:
                self.modification4 = self.dic["modification4"]
            if "enzyme" in self.dic:
                self.enzyme = self.dic["enzyme"]
            if "numMods" in self.dic:
                self.numMods = self.dic["numMods"]
            if "memory" in self.dic:
                self.memory = self.dic["memory"]
            if "mods" in self.dic:
                self.mods = self.dic["mods"]
            if "tolerance" in self.dic:
                self.tolerance = self.dic["tolerance"]
            if "fragmentation" in self.dic:
                self.fragmentation = self.dic["fragmentation"]
            if "protocol" in self.dic:
                self.protocol = self.dic["protocol"]
            if "instrument" in self.dic:
                self.instrument = self.dic["instrument"]
            if "decoy" in self.dic:
                self.decoy = self.dic["decoy"]
            if "termini" in self.dic:
                self.termini = self.dic["termini"]
            if "isotope" in self.dic:
                self.isotope = self.dic["isotope"]
            if "length" in self.dic:
                self.length = self.dic["length"]
            if "charge" in self.dic:
                self.charge = self.dic["charge"]
            if "matches" in self.dic:
                self.matches = self.dic["matches"]
            if "missCleavages" in self.dic:
                self.missCleavages = self.dic["missCleavages"]
            if "irrCleavages" in self.dic:
                self.irrCleavages = self.dic["irrCleavages"]
            if "ParentError" in self.dic:
                self.ParentError = self.dic["ParentError"]
            if "Fdr" in self.dic:
                self.Fdr = self.dic["Fdr"]

            self.get_args_fullpath('peptide_search')
            self.run_pep_search()

    def postms(self):
        """ Set input arguments for post ms mode. Calls run_postms() function to process results. """
        self.set_postms_standard()
        if "Mass_spec" not in self.dic:
            print("Please inform the folder containing your MzML files.")
        else:
            if "Mass_spec" in self.dic:
                self.Mass_spec = os.path.abspath(self.dic["Mass_spec"])
            if "Transcriptome" in self.dic:
                self.Transcriptome = self.dic["Transcriptome"]
            if "proteome" in self.dic:
                self.proteome = self.dic["proteome"]
            if "genome" in self.dic:
                self.genome = self.dic["genome"]
            if "genbank" in self.dic:
                self.genbank = self.dic["genbank"]
            self.get_args_fullpath('postms')
            self.run_postms()


    def set_wd(self):
        """ Sets working directory for the mode its applied on. """
        if not os.path.exists(os.path.abspath(self.dic["outdir"])):
            cmd_dir = "mkdir %s" % self.dic["outdir"]
            os.system(cmd_dir)
        os.chdir(self.dic["outdir"])
        print(self.dic["outdir"])

    def run_db(self):
        """ Runs the database generation step. Should be applied inside the database() function. """
        self.set_wd()
        genome = dg.OrfPrediction(self)
        genome.identify_orfs()
        print("\nORFs identified \nNow performing steps to generate the GENOME database\n")
        genome_db = dg.Database("genome_ORFs.fasta", self.proteome, "genome")
        genome_db.blast_to_Proteome()
        genome_db.extract_entries()
        cmd_cat_entries = 'cat genome_entries.txt genome_entries_short.txt > genome_entries_cat.txt'
        os.system(cmd_cat_entries)
        genome_db.remove_entries()
        genome_db.create_custom()
        if self.Transcriptome is not None:
            print("\nGenome database generated. Now performing steps to generate the TRANSCRIPTOME database.\n")
            transcriptome_db = dg.Database("transcriptome_ORFs.fasta", self.proteome, "transcriptome")
            transcriptome_db.blast_to_Proteome()
            transcriptome_db.extract_entries()
            cmd_cat_entries_rna = 'cat transcriptome_entries.txt transcriptome_entries_short.txt > ' \
                                  'transcriptome_entries_cat.txt'
            os.system(cmd_cat_entries_rna)
            transcriptome_db.remove_entries()
            transcriptome_db.create_custom()
        print("Database generation step complete. Look for databases in %s" % self.outdir)

    def run_assembly(self):
        """ Runs the assembly step. Should be applied inside the assembly() function. """
        self.set_wd()
        rna_seq = dg.Assembly(self)
        rna_seq.denovo_assembly()

    def run_pep_search(self):
        """ Runs the peptide search step. Should be called inside the peptide_search() function. """
        self.set_wd()
        genome = ps.PeptideSearch("Genome", self.Mass_spec, "genome_database.fasta", self)
        genome.peptide_identification()
        genome.peptide_filtering()
        if self.Transcriptome == "YES":
            transcriptome = ps.PeptideSearch("Transcriptome", self.Mass_spec, "transcriptome_database.fasta", self)
            transcriptome.peptide_identification()
            transcriptome.peptide_filtering()

    def run_postms(self):
        """ Runs the postMS step. Should be called inside the post_ms() function. """
        self.set_wd()
        pms.merge_databases()
        genome_results = pms.Results("Genome", "genome", "genome_database.fasta")
        genome_results.rename_files()
        genome_results.write_results()
        genome_results.merge()
        genome_results.add_orf()
        genome_results.add_orf_sequence()
        genome_results.total_spec_count()
        genome_results.genome_coordinates()
        genome_results.check_unique(self)
        genome_results.coverage()
        genome_results.summarize_results()
        genome_results.unique_orfs()
        genome_results.summarize_dna_rna()
        genome_results.create_fast("Genome")
        genome_results.minimal_results()
        genome_results.minimum_runs_compact("15", "genome")
        if self.genbank is not None:
            genome_ncrnas = pms.GenomicContext(self, "genome", "15")
            genome_ncrnas.match_features()
            genome_ncrnas.add_to_df()
        if self.Transcriptome == "YES":
            transcriptome_results = pms.Results("Transcriptome", "transcriptome",
                                                "transcriptome_database.fasta")
            transcriptome_results.rename_files()
            transcriptome_results.write_results()
            transcriptome_results.merge()
            transcriptome_results.add_orf()
            transcriptome_results.add_orf_sequence()
            transcriptome_results.total_spec_count()
            transcriptome_results.genome_coordinates()
            transcripts = pms.TranscriptLocalizer(self)
            transcripts.find_transcript()
            transcripts.fix_nucmer()
            transcripts.add_coords()
            transcriptome_results.check_unique(self)
            transcriptome_results.coverage()
            transcriptome_results.summarize_results()
            transcriptome_results.summarize_dna_rna()
            pms.create_venn()
            transcriptome_results.add_names()
            genome_results.create_fast("Both")
            transcriptome_results.create_fast("Rna")
            transcriptome_results.merge_both_results()
            both_results = pms.Results("Genome", "Both", "_")
            both_results.minimal_results()
            both_results.minimum_runs_compact("15", "both")
            transcriptome_results.minimal_results()
            transcriptome_results.minimum_runs_compact("15", "transcriptome")
            if self.genbank is not None:
                rna_features = pms.GenomicContext(self, "transcriptome", "15")
                rna_features.match_features()
                rna_features.add_to_df()
                both_features = pms.GenomicContext(self, "both", "15")
                both_features.match_features()
                both_features.add_to_df()
