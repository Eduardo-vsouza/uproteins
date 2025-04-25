import os
import pathlib
import subprocess

import cli


class Error(Exception):
    pass


class ReadError(Error):
    def __init__(self, message="Inconsistent read information. Provide a set of paired-end reads, or a set of single"
                               " end reads."):
        self.message = message
        super().__init__(self.message)


class ReadMapper(object):
    def __init__(self, args):
        self.args = args
        self.gtf = args.gtf
        self.genome = args.genome
        self.is_paired = self._check_read_type()
        self.strandness = self.args.strandness

    def _check_read_type(self):
        # We assume validate_assembly made sure that only --single OR
        # (--reads1 --reads2) was given
        if self.args.single is None:
            return True
        return False

    def create_indexes(self):
        cmd_index = f'hisat2-build {self.genome} genome'
        os.system(cmd_index)
        return self

    def align_reads(self):
        pathlib.Path('HISAT/sam').mkdir(exist_ok=True, parents=True)

        reads = self.args.single
        if self.is_paired:
            reads = zip(self.args.reads1, self.args.reads2)

        # read_data is either a tuple of paths (-1 and -2) or a single path
        for i, read_data in enumerate(reads):
            read_name = f'aligned_{i}'

            hisat_strandness = f'--rna-strandness {self.strandness} ' \
                if self.strandness is not None else ''

            hisat_input = f'-U {read_data}'
            if self.is_paired:
                read_one, read_two = read_data
                hisat_input = f'-1 {read_one} -2 {read_two}'

            cmd_hisat = (
                f'hisat2 -x genome {hisat_input} --dta {hisat_strandness} '
                f'-S HISAT/sam/{read_name}.sam'
            )

            result = subprocess.run(cmd_hisat)
            if result.returncode != 0:
                cli.exit(
                    3,
                    'hisat2 finished with non-zero return value: '
                    f'{result.returncode}'
                )

            self._sort_alignment(read_name)

        return self

    def _sort_alignment(self, read_name: str):
        """ Uses samtools to sort the alignments in the SAM file. """
        pathlib.Path("HISAT/bam").mkdir(exist_ok=True, parents=True)
        cmd_sam = f'samtools view -Su HISAT/sam/{read_name}.sam | samtools sort -o HISAT/bam/{read_name}.bam'
        os.system(cmd_sam)
        cmd_rv = f'rm HISAT/sam/{read_name}.sam'
        os.system(cmd_rv)
        return self


class TranscriptAssembly(object):
    def __init__(self, args):
        self.args = args
        self.threads = args.threads
        self.gtf = self.args.gtf
        self.strandness = self._check_strand()

    def _check_strand(self):
        if self.args.strandness == "RF":
            # or `self.args.strandness == 'R'`?
            strandness = "--rf "
        elif self.args.strandness == "FR" or self.args.strandness == 'F':
            strandness = "--fr "
        else:
            strandness = ''
        return strandness

    def assemble(self):
        """ Runs StringTie to assemble transcripts """
        files = os.listdir("HISAT/bam")
        if not os.path.exists("StringTie"):
            cmd_dir = 'mkdir StringTie'
            os.system(cmd_dir)
        if not os.path.exists("StringTie/gtf_intermediates"):
            cmd_dir2 = 'mkdir StringTie/gtf_intermediates'
            os.system(cmd_dir2)
        for file in files:
            # self.strandness already has a trailing space
            cmd_str = f'stringtie HISAT/bam/{file} -o StringTie/gtf_intermediates/{file[:-3]}gtf -p {self.threads} -G {self.gtf} ' \
                      f'{self.strandness}-l uproteins -m 30 -A StringTie/gene_abundances.txt'
            os.system(cmd_str)
        return self

    def create_gtf_list(self):
        files = os.listdir("StringTie/gtf_intermediates")
        all_gtfs = []
        for file in files:
            if '.gtf' in file:
                all_gtfs.append("StringTie/gtf_intermediates/"+file+"\n")
        with open("StringTie/gtf_list.txt", 'w') as out:
            out.writelines(all_gtfs)
        return self

    def merge_transcripts(self):
        cmd_str = f'stringtie --merge StringTie/gtf_list.txt -o assembled.gtf -G {self.gtf} -m 30 -l uproteins'
        os.system(cmd_str)
        gtf_path = 'assembled.gtf'
        return gtf_path


class CompareTranscripts(object):
    def __init__(self, args, assembled_gtf):
        self.assembled_gtf = assembled_gtf
        self.args = args
        self.__check_dir()
        self.path = self.__check_path()
        self.gffreadPath = self.__check_gffread_path()
        self.copy_genome()
        self.index_genome()

    @staticmethod
    def __check_dir():
        if not os.path.exists("RNA_comparisons"):
            cmd_dir = 'mkdir RNA_comparisons'
            os.system(cmd_dir)

    def copy_genome(self):
        cmd_cp = f'cp {self.args.genome} ./genome.fasta'
        os.system(cmd_cp)

    def __check_path(self):
        if self.args.gff_compare_path is not None:
            gff_path = self.args.gff_compare_path
        else:
            gff_path = 'gffcompare'
        return gff_path

    def index_genome(self):
        cmd_sam = f'samtools faidx genome.fasta'
        os.system(cmd_sam)

    def run_gffcompare(self):
        cmd_gff = f'{self.path} -r {self.args.gtf} {self.assembled_gtf} -o RNA_comparisons/compared_assembly -s' \
                  f' genome.fasta'
        os.system(cmd_gff)
        return 'RNA_comparisons'

    def __check_gffread_path(self):
        if self.args.gffread_path is not None:
            gffread = self.args.gffread_path
        else:
            gffread = 'gffread'
        return gffread

    def extract_sequences(self):
        cmd_read = f'{self.gffreadPath} -w HISAT/transcripts.fasta -g genome.fasta {self.assembled_gtf}'
        os.system(cmd_read)
        return 'HISAT/transcripts.fasta'


# class Transcript(object):
#     def __init__(self, seqname, start, end, strand, attributes):
#         self.name = seqname
#         self.start = start
#         self.end = end
#         self.strand = strand
#         self.attrs = attributes


# class AssemblyExtractor(object):
#     def __init__(self, assembled_gtf):
#         # self.args = args
#         self.genome = 'home/eduardo/Documents/smorfs_smeg_stringtie/rna-seq/annotation_files/genome.fasta'
#         self.gtf = assembled_gtf
#         self.df = self.__read_gtf()
#         self.transcripts = []
#
#     def __read_gtf(self):
#         df = pd.read_csv(self.gtf, sep="\t", header=2)
#         df.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
#         transcript_df = df[df["feature"] == "transcript"]
#         return transcript_df
#
#     def get_transcripts(self):
#         seqnames = self.df["seqname"].tolist()
#         starts = self.df["start"].tolist()
#         ends = self.df["end"].tolist()
#         strands = self.df["strand"].tolist()
#         attrs = self.df["attributes"].tolist()
#         for i in range(len(seqnames)):
#             transcript = Transcript(seqnames[i], int(starts[i])-1, int(ends[i])-1, strands[i], attrs[i])
#             self.transcripts.append(transcript)
#         return self
#
#     def lookup_transcripts(self):
#         records = SeqIO.parse(self.genome, 'fasta')
#         fasta = []
#         for rna in self.transcripts:
#             for record in records:
#                 if rna.seqname in record.id:
#                     seq = record.seq[rna.start:rna.end]
#                     gene_id = rna.attrs.split(";")
#                     print(gene_id)
#                     break
