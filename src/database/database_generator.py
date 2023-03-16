import os
import sys

from Bio import SeqIO
from Bio.Blast import NCBIXML

from ..translate import GenomeReader as tr


path = sys.path[0]


class Assembly(object):
    def __init__(self, args):
        self.args = args
        """ Performs DE NOVO Transcriptome Assembly using Trinity if <-Transcriptome YES> is specified in the command
            line. """

    def denovo_assembly(self):
        reads1 = "--left"
        reads2 = "--right"
        single_reads = "--single"
        lib = ""
        cpu = ""
        memory = ""
        trinity_path = f"{path}/dependencies/trinity_for_uproteins"
        if self.args.lib_type is not None:
            lib = " --SS_lib_type " + self.args.lib_type
        if self.args.cpu is not None:
            cpu = " --CPU {}".format(self.args.cpu)
        if self.args.memory is not None:
            memory = " --max_memory {}".format(self.args.memory)
        if self.args.rna_right is not None:
            Cmd_T_0 = f'{trinity_path}/Trinity --no_salmon --seqType fq %s %s %s %s ' \
                      '--min_contig_length 30%s%s%s --output trinity/' % (reads1, self.args.rna_left, reads2,
                                                                       self.args.rna_right, lib, cpu, memory)
            os.system(Cmd_T_0)
        if self.args.rna_single is not None:
            Cmd_T_1 = f'{trinity_path}/Trinity --seqType fq %s %s --min_contig_' \
                      'length 30%s%s%s --output trinity/' % (single_reads, self.args.rna_single, lib, cpu, memory)
            os.system(Cmd_T_1)

    def genome_mapping(self):
        """ This is deprecated. """
        build = 'gmap_build -D ./ -d %s trinity/Trinity.fasta' % self.args.genome
        os.system(build)
        gmap = 'gmap -Z -d %s trinity/Trinity.fasta > mapped_assembly' % self.args.genome
        os.system(gmap)


class OrfPrediction(object):
    def __init__(self, args):
        self.args = args
        self.fix_genome_entry()

    def fix_genome_entry(self):
        genome = self.args.genome
        fix_genome = []
        records = SeqIO.parse(genome, 'fasta')
        i = 0
        for record in records:
            if i > 0:
                fix_genome.append(f">NC_{self.args.organism}_{i}\n{record.seq}\n")
            else:
                fix_genome.append(f">NC_{self.args.organism}\n{record.seq}\n")
            i += 1
        with open('fixed_genome.fasta', 'w') as fa:
            fa.writelines(fix_genome)
        self.args.genome = os.path.abspath('fixed_genome.fasta')

    def identify_orfs(self):
        genome = tr(self.args.genome, starts=self.args.starts, ends=self.args.ends, filetype="genome", args=self.args)
        genome.get_sequence()
        genome.three_frame_tr(genome.entries, genome.sequence, "+")
        genome.three_frame_tr(genome.c_entries, genome.c_seqs, "-")
        genome.translate()
        if self.args.Transcriptome is not None:
            transcriptome = tr("trinity/Trinity.fasta", starts=self.args.starts, ends=self.args.ends,
                                            filetype="transcriptome", args=self.args)
            transcriptome.get_sequence()
            transcriptome.three_frame_tr(transcriptome.entries, transcriptome.sequence, "+")
            transcriptome.translate()

    # def Identify_ORFs(self):
    #     """ Translates the genome and transcriptome to the 6 and 3 reading frames, respectively. """
    #     tab = 0
    #     circ = "NO"
    #     max_size = "1000000"
    #     if self.args.maxsize is not None:
    #         max_size = self.args.maxsize
    #     if self.args.table is not None:
    #         tab = self.args.table
    #     if self.args.circular is not None:
    #         circ = self.args.circular
    #     if self.args.Transcriptome is not None:
    #         Cmd_T_2 = 'getorf -sequence trinity/Trinity.fasta -outseq primary_RNA_ORFs.fasta -table %s -scircular1 NO ' \
    #                   '-find 1 -methionine NO -reverse NO -minsize 10 -maxsize %s' % (tab, max_size)
    #         os.system(Cmd_T_2)
    #         with open("transcriptome_ORFs.fasta", 'w') as replaced, open("primary_RNA_ORFs.fasta", 'r') as original:
    #             fasta = original.readlines()
    #             new = ""
    #             for line in fasta:
    #                 new = new + line.replace(" ", "_")
    #             replaced.writelines(new + "\n")
    #     Cmd0 = 'getorf -sequence %s -outseq primary_ORFs.fasta -table %s -scircular1 %s -find 1 ' \
    #            '-methionine NO -minsize 10 -maxsize %s' % (self.args.genome, tab, circ, max_size)
    #     os.system(Cmd0)
    #     with open("genome_ORFs.fasta", 'w') as replaced, open("primary_ORFs.fasta", 'r') as original:
    #         fasta = original.readlines()
    #         new = ""
    #         for line in fasta:
    #             new = new + line.replace(" ", "_")
    #         replaced.writelines(new + "\n")


class Database(object):
    def __init__(self, ob, p, filetype):
        self.orf_to_blast = ob
        self.proteome = p
        self.filetype = filetype
        self.blast_dir = f'{path}/dependencies/blast_for_uproteins/bin'

    def mark_annotated(self):
        # if not os.path.exists('annotated.fasta'):
        marked = []
        records = SeqIO.parse(self.proteome, 'fasta')
        for record in records:
            marked.append(f'>{str(record.description)}_ANNO\n{str(record.seq)}\n')
        with open('annotated.fasta', 'w') as out:
            out.writelines(marked)
        self.annotated = 'annotated.fasta'

    def unify(self):
        cmd_cat = f'cat {self.orf_to_blast} {self.annotated} > {self.filetype}_database.fasta'
        os.system(cmd_cat)


    def blast_to_Proteome(self):
        """ Aligns the ORFs to the annotated proteome with Blastp in order to identify annotated entries. """
        Cmd_short = f'blastp -query %s -subject %s -outfmt 5 ' \
                    '-qcov_hsp_perc 90 -comp_based_stats 0 -task blastp-short -out %s_Blasted_to_Uniprot_short.xml' \
                    % (self.orf_to_blast, self.proteome, self.filetype)
        os.system(Cmd_short)
        cmd_blast = f'blastp -query %s -subject %s -outfmt 5 -qcov_hsp_perc 90 ' \
                    '-out %s_Blasted_to_Uniprot.xml' % (self.orf_to_blast, self.proteome, self.filetype)
        os.system(cmd_blast)

    def extract_entries(self):
        """ Writes the annotated entries into another file. """
        genome_handle = open('%s_Blasted_to_Uniprot.xml' % self.filetype, 'r')
        blast_records_genome = NCBIXML.parse(genome_handle)
        with open("%s_entries.txt" % self.filetype, "w") as resss:
            results = []
            try:

                for record in blast_records_genome:
                    for alignment in record.alignments:
                        for hsp in alignment.hsps:
                            if hsp.identities >= (hsp.align_length - 2):
                                results.append(record.query + "\n")
            except:
                print("foo")
            resss.writelines(results)
        genome_handle_short = open('%s_Blasted_to_Uniprot_short.xml' % self.filetype, 'r')
        blast_records_genome_short = NCBIXML.parse(genome_handle_short)
        with open("%s_entries_short.txt" % self.filetype, "w") as resss_short:
            results = []
            try:
                for record in blast_records_genome_short:
                    for alignment in record.alignments:
                        for hsp in alignment.hsps:
                            if hsp.identities >= (hsp.align_length - 2):
                                results.append(record.query + "\n")
            except:
                print('foo')
            resss_short.writelines(results)
        cmd_cat_entries = 'cat %s_entries.txt %s_entries_short.txt > %s_entries_cat.txt' % (self.filetype, self.filetype,
                                                                                            self.filetype)
        os.system(cmd_cat_entries)

    def remove_entries(self):
        """ Removes the annotated ORFs identified in the previous blastp step from the ORFs database. """
        entries = set([])
        lines = open("%s_entries_cat.txt" % self.filetype, 'r')
        for l in lines:
            l = l.strip()
            entries.add(str(l).replace('>', ""))
        with open("%s" % self.orf_to_blast, 'r') as original, open('%s_database_no_anno.fasta' % self.filetype, 'w') \
                as removed:
            records = SeqIO.parse(original, 'fasta')
            print("Removing database entries that are already annotated.")
            for record in records:
                if record.id not in entries:
                    SeqIO.write(record, removed, 'fasta')

    def create_custom(self):
        """ Merges both predicted DB (composed by unannotated smORFs) and RefSeq DB into a single, custom database. """
        cmd_cat_db = 'cat %s_database_no_anno.fasta %s > %s_database.fasta' % (self.filetype, self.proteome,
                                                                             self.filetype)
        os.system(cmd_cat_db)
