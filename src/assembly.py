import os

from Bio import SeqIO
import pandas as pd


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
        self.__check_read_type()
        self.libtype = None
        self.strandness = self.args.strandness

    def __check_read_type(self):
        if self.args.reads1 is not None and self.args.single is None:
            self.__get_paired_reads()
        elif self.args.reads1 is None and self.args.single is not None:
            self.__get_single_reads()
        # elif self.args.reads1 is None and self.args.singles is None:
        #     raise ReadError
        else:
            raise ReadError

    def __get_single_reads(self):
        reads = self.args.single.split(",")
        read_list = []
        # for read in reads:
        #         read_list.append(os.path.abspath(read))
        self.libtype = "single"
        self.reads = reads

    def __get_paired_reads(self):
        reads1 = self.args.reads1.split(",")
        reads2 = self.args.reads2.split(",")
        self.libtype = "paired"
        self.reads1 = reads1
        self.reads2 = reads2

    def create_indexes(self):
        cmd_index = f'hisat2-build {self.genome} genome'
        os.system(cmd_index)
        return self

    def align_reads(self):
        if not os.path.exists("HISAT"):
            cmd_dir = 'mkdir HISAT'
            os.system(cmd_dir)
        if not os.path.exists("HISAT/sam"):
            cmd_dir2 = 'mkdir HISAT/sam'
            os.system(cmd_dir2)
        # if self.libtype == "single":
        print('alright')
        i = 0
        for read in self.reads:
            read_name = f'aligned_{i}'
            cmd_hisat = f'hisat2 -x genome -U {read} --dta --rna-strandness {self.strandness} -S HISAT/sam/{read_name}.sam'
            os.system(cmd_hisat)
            self.__sort_aln(i)
            i += 1
        # elif self.libtype == "paired":
        #     i = 0
        #     for read1, read2 in zip(self.reads1, self.reads2):
        #         i += 1
        #         name1 = read1.find("_")
        #         read1_name = f"aligned_{i}"
        #         cmd_hisat = f'hisat2 -x genome -1 {read1} -2 {read2} --dta --rna-strandness {self.strandness} -S' \
        #                     f' HISAT/sam/{read1_name}.sam'
        #         os.system(cmd_hisat)
        #         self.sort_aln(i)
        return self

    def __sort_aln(self, i):
        """ Uses samtools to sort the alignments in the SAM file. """
        # files = os.listdir("HISAT/sam")
        if not os.path.exists("HISAT/bam"):
            cmd_dir = 'mkdir HISAT/bam'
            os.system(cmd_dir)

        # for file in files:
        cmd_sam = f'samtools view -Su HISAT/sam/aligned_{i}.sam | samtools sort -o HISAT/bam/aligned_{i}.bam'
        os.system(cmd_sam)
        cmd_rv = f'rm HISAT/sam/aligned_{i}.sam'
        os.system(cmd_rv)
        return self


class TranscriptAssembly(object):
    def __init__(self, args):
        self.args = args
        self.threads = args.threads
        self.gtf = self.args.gtf
        self.strandness = self.check_strand()

    def check_strand(self):
        if self.args.strandness == "RF":
            strandness = "--rf "
        elif self.args.strandness == "FR":
            strandness = "--fr "
        elif self.args.strandness == 'F':
            strandness = "--fr "
        else:
            strandness = self.args.strandness
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
            gff_path = '/data/gffcompare-0.12.6.Linux_x86_64/gffcompare'
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

#
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
