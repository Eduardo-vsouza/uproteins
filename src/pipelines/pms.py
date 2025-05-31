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


from ..postprocess import PostPercolator, ExtendedInformation, PercolatorProcessing, AllSub, TSVConverter, ResultsWrapper
from ..sequtils.orflib import AltCodons
from ..upstream import SDInspection
from ..sequtils.__helpers import FiletypeError


class PostMSPipeline(object):
    def __init__(self, args, filetype, folder, qvalue=0.01, testing=False):
        self.args = args
        self.filetype = filetype
        self.folder = folder
        self.qValue = qvalue
        self.testing = testing

    def run(self):
        self._run_percolator()
        self._process_percolator()
        self._select_codons()
        self._reformat_results()

    def _run_percolator(self):
        perc = PercolatorProcessing(folder=self.folder, filetype=self.filetype)
        perc.create_metafiles().convert_to_pin()
        perc.percolate()

    def _process_percolator(self):
        tsv = TSVConverter(self.folder)
        tsv.convert_files()

        data_filter = PostPercolator(self.args, folder=self.folder, filetype=self.filetype)
        data_filter.convert_output()
        if self.filetype == 'genome':
            data_filter.get_coordinates_genome()
        elif self.filetype == 'transcriptome':
            data_filter.get_coordinates_rna()
        else:
            raise FiletypeError
        data_filter.filter_novel()
        data_filter.unique_peptides()
        data_filter.msgf_info()
        data_filter.protein_seqs()
        data_filter.add_coordinates(qvalue=self.qValue)

    def _select_codons(self):
        if self.filetype == 'genome':
            alts_pre_rf = AltCodons(file='Genome/post_perc/genome_results_02.txt', genome=self.args.genome,
                                    maxsize=self.args.maxsize, testing=self.testing)
        elif self.filetype == 'transcriptome':
            alts_pre_rf = AltCodons(file='Transcriptome/post_perc/transcriptome_results_02.txt',
                                    genome=self.args.genome, maxsize=self.args.maxsize,
                                    transcriptome_gff='assembled.gtf', assembly="HISAT/transcripts.fasta",
                                    subset="Transcriptome", testing=self.testing)
        else:
            raise FiletypeError
        alts_pre_rf.extend_orfs(args=self.args)
        if self.args.rrna is not None:
            if self.filetype == 'genome':
                genome_rbs = SDInspection(self.args, filetype="genome", folder="Genome",
                                          alternatives=alts_pre_rf.alternatives, testing=self.testing)
            elif self.filetype == 'transcriptome':
                genome_rbs = SDInspection(self.args, filetype='transcriptome', folder='Transcriptome',
                                          alternatives=alts_pre_rf.alternatives,
                                          transcriptome_gff='assembled.gtf',
                                          transcripts="HISAT/transcripts.fasta", subset="Transcriptome",
                                          testing=self.testing)
            else:
                raise FiletypeError
            alts = genome_rbs.get_free_energy()
            alts_pre_rf.alternatives = alts
        alts_pre_rf.sort_by_coordinates()
        alts_pre_rf.sort_by_atg()
        alts_pre_rf.sort_by_shine()
        alts_pre_rf.sort_by_peptides()
        priorities = alts_pre_rf.get_priorities()
        ext = ExtendedInformation(folder=self.folder, filetype=self.filetype, alternatives=alts_pre_rf.alternatives)
        ext.filter_alternatives(priorities)
        ext.extract_spectra()

    def _reformat_results(self):
        results = ResultsWrapper(df=f'{self.folder}/post_perc/{self.filetype}_results_04.txt', folder=self.folder,
                                 filetype=self.filetype)
        results.reformat(pre_validation=True)

