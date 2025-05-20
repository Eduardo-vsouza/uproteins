import pathlib
from importlib import resources as rsrc

from src import uproteins
from tests import resources


def test_full_run(tmp_path):
    assembly_args = [
        '--outdir', str(tmp_path),
        'assembly',
        '--single',
            'tests/resources/ERR262980.fastq',
            'tests/resources/ERR262982.fastq',
            'tests/resources/ERR262983.fastq',
        '--genome', 'tests/resources/genome.fasta',
        '--gtf', 'tests/resources/mtb.gtf',
        '--strandness', 'F',
    ]
    database_args = [
        '--outdir', str(tmp_path),
        'database',
        '--genome', 'tests/resources/genome.fasta',
        '--proteome', 'tests/resources/proteome.fasta',
        '--starts', 'ATG', 'GTG', 'TTG', 'CTG',
        '--minsize', '30',
        '--maxsize', '300',
        '--transcriptome'
    ]
    ms_args = [
        '--outdir', str(tmp_path),
        'ms',
        '--mass-spec', 'tests/resources/mzml/',
        '--inst', '2',
        '--t', '800ppm',
        '--transcriptome'
    ]
    postms_args = [
        '--outdir', str(tmp_path),
        'postms',
        '--genome', 'tests/resources/genome.fasta',
        '--proteome', 'tests/resources/proteome.fasta',
        '--rrna', 'tests/resources/rrna.fna',
        '--gff', 'tests/resources/mtb.gtf',
        '--mass-spec', 'tests/resources/mzml/',
        '--transcriptome'
    ]
    validate_args = [
        '--outdir', str(tmp_path),
        'validate',
        '--transcriptome'
    ]
    uproteins(assembly_args)
    uproteins(database_args)
    uproteins(ms_args)
    uproteins(postms_args)
    uproteins(validate_args)

    # Paths to the output
    genome_pre: pathlib.Path = (
        tmp_path / 'Genome' / 'genome_pre_validation_results.txt'
    )
    genome_post: pathlib.Path = (
        tmp_path / 'Genome' / 'genome_post_validation_results.txt'
    )
    
    transcriptome_pre: pathlib.Path = (
        tmp_path
        / 'Transcriptome'
        / 'transcriptome_pre_validation_results.txt'
    )
    transcriptome_post: pathlib.Path = (
        tmp_path
        / 'Transcriptome'
        / 'transcriptome_post_validation_results.txt'
    )

    # Make sure the result files were created
    assert genome_pre.is_file()
    assert genome_post.is_file()
    assert transcriptome_pre.is_file()
    assert transcriptome_post.is_file()

    # And that they have the expect results
    with rsrc.path(resources, 'results') as results:
        ok_genome_pre = (
            results
            / 'Genome'
            / 'genome_pre_validation_results.txt'
        )
        ok_genome_post = (
            results
            / 'Genome'
            / 'genome_post_validation_results.txt'
        )
        
        ok_transcriptome_pre = (
            results
            / 'Transcriptome'
            / 'transcriptome_pre_validation_results.txt'
        )
        ok_transcriptome_post = (
            results
            / 'Transcriptome'
            / 'transcriptome_post_validation_results.txt'
        )

        assert genome_pre.read_text() == ok_genome_pre.read_text()
        assert (
            transcriptome_pre.read_text()
            == ok_transcriptome_pre.read_text()
        )

        assert genome_post.read_text() == ok_genome_post.read_text()
        assert (
            transcriptome_post.read_text()
            == ok_transcriptome_post.read_text()
        )
