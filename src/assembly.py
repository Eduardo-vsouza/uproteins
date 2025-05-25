import os
import pathlib
import shutil
import subprocess

from src import cli


class Error(Exception):
    pass


class ReadError(Error):
    def __init__(
        self,
        message="Inconsistent read information. Provide a set of paired-end "
        "reads, or a set of single end reads."
    ):
        self.message = message
        super().__init__(self.message)


class ReadMapper:
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
        cmd_index = f'hisat2-build {self.genome} genome'.split(' ')
        result = subprocess.run(cmd_index)
        if (code := result.returncode) != 0:
            err = f'hisat2-build finished with non-zero return value: {code}'
            cli.exit(3, err)
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
            ).split(' ')

            result = subprocess.run(cmd_hisat)
            if (code := result.returncode) != 0:
                err = f'hisat2 finished with non-zero return value: {code}'
                cli.exit(3, err)

            self._sort_alignment(read_name)

        return self

    def _sort_alignment(self, read_name: str):
        """Use samtools to sort the alignments in the SAM file."""
        pathlib.Path("HISAT/bam").mkdir(exist_ok=True, parents=True)

        # run samtools view and pipe it to samtools sort
        cmd_view = f'samtools view -Su HISAT/sam/{read_name}.sam'.split(' ')
        view_result = subprocess.run(
            cmd_view,
            stdout=subprocess.PIPE,
        )
        if (code := view_result.returncode) != 0:
            err = f'samtools view finished with non-zero return value: {code}'
            cli.exit(3, err)

        cmd_sort = f'samtools sort -o HISAT/bam/{read_name}.bam'.split(' ')
        sort_result = subprocess.run(
            cmd_sort,
            input=view_result.stdout,
        )
        if (code := sort_result.returncode) != 0:
            err = f'samtools sort finished with non-zero return value: {code}'
            cli.exit(3, err)

        # cleanup
        pathlib.Path(f'HISAT/sam/{read_name}.sam').unlink()
        return self


class TranscriptAssembly:
    def __init__(self, args):
        self.args = args
        self.threads = args.threads
        self.gtf = self.args.gtf
        self.strandness = self._check_strand()

    def _check_strand(self):
        if self.args.strandness == "RF":
            # or `self.args.strandness == 'R'`?
            return "--rf "
        elif self.args.strandness == "FR" or self.args.strandness == 'F':
            return "--fr "
        return ''

    def assemble(self):
        """Run StringTie to assemble transcripts."""
        path = pathlib.Path('StringTie/gtf_intermediates')
        path.mkdir(exist_ok=True, parents=True)
        files = os.listdir("HISAT/bam")
        for file in files:
            # self.strandness already has a trailing space
            cmd_str = (
                f'stringtie HISAT/bam/{file} -o '
                f'StringTie/gtf_intermediates/{file[:-3]}gtf -p '
                f'{self.threads} -G {self.gtf} {self.strandness}-l uproteins '
                '-m 30 -A StringTie/gene_abundances.txt'
            ).split(' ')
            result = subprocess.run(cmd_str)
            if (code := result.returncode) != 0:
                err = f'stringtie finished with non-zero return value: {code}'
                cli.exit(3, err)
        return self

    def create_gtf_list(self):
        files = os.listdir("StringTie/gtf_intermediates")
        all_gtfs = []
        for file in files:
            if '.gtf' in file:
                all_gtfs.append(f"StringTie/gtf_intermediates/{file}\n")
        with open("StringTie/gtf_list.txt", 'w') as out:
            out.writelines(all_gtfs)
        return self

    def merge_transcripts(self):
        cmd_str = (
            'stringtie --merge StringTie/gtf_list.txt -o assembled.gtf -G '
            f'{self.gtf} -m 30 -l uproteins'
        ).split(' ')
        result = subprocess.run(cmd_str)
        if (code := result.returncode) != 0:
            err = f'stringtie finished with non-zero return value: {code}'
            cli.exit(3, err)
        gtf_path = 'assembled.gtf'
        return gtf_path


class CompareTranscripts:
    def __init__(self, args, assembled_gtf):
        self.assembled_gtf = assembled_gtf
        self.args = args
        self.gffcompare_path = self.args.gffcompare_path
        self.gffread_path = self.args.gffread_path
        self.copy_genome()
        self.index_genome()

    def copy_genome(self) -> None:
        shutil.copy(self.args.genome, './genome.fasta')

    def index_genome(self):
        cmd_sam = 'samtools faidx genome.fasta'.split(' ')
        result = subprocess.run(cmd_sam)
        if (code := result.returncode) != 0:
            err = f'samtools faidx finished with non-zero status code: {code}'
            cli.exit(3, err)

    def run_gffcompare(self):
        pathlib.Path('RNA_comparisons').mkdir(exist_ok=True)
        cmd_compare = (
            f'{self.gffcompare_path} -r {self.args.gtf} '
            f'{self.assembled_gtf} -o '
            'RNA_comparisons/compared_assembly -s genome.fasta'
        ).split(' ')
        result = subprocess.run(cmd_compare)
        if (code := result.returncode) != 0:
            err = f'gffcompare finished with non-zero return code: {code}'
            cli.exit(3, err)
        return 'RNA_comparisons'

    def extract_sequences(self):
        cmd_read = (
            f'{self.gffread_path} -w HISAT/transcripts.fasta -g '
            f'genome.fasta {self.assembled_gtf}'
        ).split(' ')
        result = subprocess.run(cmd_read)
        if (code := result.returncode) != 0:
            err = f'gffread finished with non-zero return code: {code}'
            cli.exit(3, err)
        return 'HISAT/transcripts.fasta'
