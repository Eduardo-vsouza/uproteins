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


class MSArgs(object):
    def __init__(self, outdir):
        self.modification1 = "IAM,57.021460,C,opt,any"
        self.mods = 4
        self.enzyme = 1
        self.tolerance = 50
        self.fragmentation = 1
        self.protocol = 0
        self.instrument = 0
        self.decoy = "TRUE"
        self.termini = 1
        self.isotope = "0,1"
        self.length = "6,40"
        self.charge = '2,3'
        self.matches = 1
        self.numMods = 3
        self.missCleavages = 2
        self.irrCleavages = 1
        self.ParentError = 30
        self.Fdr = 0.01
        self.outdir = os.path.abspath(outdir)
        self.Transcriptome = "YES"
        self.memory = "6000"
        # self.tasks = 3

    def get_ms_folder(self, test_path):
        Mass_spec = f"{test_path}/mzML"
        return Mass_spec


class TestArgs(object):
    def __init__(self, test_path, outdir, skip_assembly, skip_db, skip_ms, skip_postms):
        self.path = test_path
        self.Transcriptome = "YES"
        self.genome = f"{self.path}/genome.fasta"
        self.proteome = f"{self.path}/proteome.fasta"
        self.genbank = f"{self.path}/genbank.gb"
        self.starts = "ATG,GTG,TTG"
        self.ends = "TAA,TAG,TGA"
        self.Mass_spec = f"{self.path}/mzML"
        self.outdir = os.path.abspath(outdir)

        ## testing
        self.skip_assembly = skip_assembly
        self.skip_db = skip_db
        self.skip_ms = skip_ms
        self.skip_postms = skip_postms

        ## assembly ##
        # srr8373230
        self.rna_left = f"{self.path}/reads_1.fastq"
        self.rna_right = f"{self.path}/reads_2.fastq"
        self.lib_type = "FR"
        self.cpu = None
        self.memory = "6G"
        self.rna_single = None

        ## database ##
        self.organism = "Test_bacteria"

        ## peptide search ##
        self.msArgs = MSArgs(outdir)