import pathlib
from importlib import resources as rsrc

from src import uproteins
from tests import resources


# Since the validate method requires a venv, it is not being run at the moment
def test_full_run(tmp_path):
    genome = rsrc.files(resources).joinpath("genome.fasta")
    proteome = rsrc.files(resources).joinpath("proteome.fasta")
    gtf = rsrc.files(resources).joinpath("mtb.gtf")
    mzml = rsrc.files(resources).joinpath("mzml")
    read1 = rsrc.files(resources).joinpath("ERR262980.fastq")
    read2 = rsrc.files(resources).joinpath("ERR262982.fastq")
    read3 = rsrc.files(resources).joinpath("ERR262983.fastq")
    rrna = rsrc.files(resources).joinpath("rrna.fna")

    assembly_args = [
        '--outdir', str(tmp_path),
        'assembly',
        '--single',
            str(read1),
            str(read2),
            str(read3),
        '--genome', str(genome),
        '--gtf', str(gtf),
        '--strandness', 'F',
    ]
    uproteins(assembly_args)

    database_args = [
        '--outdir', str(tmp_path),
        'database',
        '--genome', str(genome),
        '--proteome', str(proteome),
        '--starts', 'ATG', 'GTG', 'TTG', 'CTG',
        '--minsize', '30',
        '--maxsize', '300',
        '--transcriptome'
    ]
    uproteins(database_args)

    ms_args = [
        '--outdir', str(tmp_path),
        'ms',
        '--mass-spec', str(mzml),
        '--inst', '2',
        '--t', '800ppm',
        '--transcriptome'
    ]
    uproteins(ms_args)


    postms_args = [
        '--outdir', str(tmp_path),
        'postms',
        '--genome', str(genome),
        '--proteome', str(proteome),
        '--rrna', str(rrna),
        '--gff', str(gtf),
        '--mass-spec', str(mzml),
        '--transcriptome'
    ]
    uproteins(postms_args)

    validate_args = [
        '--outdir', str(tmp_path),
        'validate',
        '--transcriptome'
    ]
    # uproteins(validate_args)

    # Paths to the output
    genome_pre: pathlib.Path = (
        tmp_path
        / 'Genome'
        / 'Results'
        / 'genome_pre_validation_results.txt'
    )
    genome_post: pathlib.Path = (
        tmp_path
        / 'Genome'
        / 'Results'
        / 'genome_post_validation_results.txt'
    )
    
    transcriptome_pre: pathlib.Path = (
        tmp_path
        / 'Transcriptome'
        / 'Results'
        / 'transcriptome_pre_validation_results.txt'
    )
    transcriptome_post: pathlib.Path = (
        tmp_path
        / 'Transcriptome'
        / 'Results'
        / 'transcriptome_post_validation_results.txt'
    )

    # Make sure the result files were created
    assert genome_pre.is_file()
    # assert genome_post.is_file()
    assert transcriptome_pre.is_file()
    # assert transcriptome_post.is_file()

    # And that they have the expect results
    with rsrc.path(resources, 'results') as results:
        ok_genome_pre = (
            results
            / 'Genome'
            / 'Results'
            / 'genome_pre_validation_results.txt'
        )
        ok_genome_post = (
            results
            / 'Genome'
            / 'Results'
            / 'genome_post_validation_results.txt'
        )
        
        ok_transcriptome_pre = (
            results
            / 'Transcriptome'
            / 'Results'
            / 'transcriptome_pre_validation_results.txt'
        )
        ok_transcriptome_post = (
            results
            / 'Transcriptome'
            / 'Results'
            / 'transcriptome_post_validation_results.txt'
        )

        assert genome_pre.read_text() == ok_genome_pre.read_text()
        assert (
            transcriptome_pre.read_text()
            == ok_transcriptome_pre.read_text()
        )

        # assert genome_post.read_text() == ok_genome_post.read_text()
        # assert (
        #     transcriptome_post.read_text()
        #     == ok_transcriptome_post.read_text()
        # )
