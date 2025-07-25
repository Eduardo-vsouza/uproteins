# Copyright © 2025 Eduardo Vieira de Souza
# Copyright © 2025 Adriana Canedo
# Copyright © 2025 Cristiano Valim Bizarro
# Copyright © 2025 Bruno Maestri A Becker
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

from argparse import Namespace
import pathlib
from importlib import resources as rsrc
import shutil

import pandas as pd
import pytest

from src.postprocess.specfilt import PostPercolator
from src import cli
from tests import resources


@pytest.fixture
def postms_dir(tmp_path_factory: pytest.TempPathFactory) -> pathlib.Path:
    """Construct and return the path to a temporary directory in the format:
    .. code-block:: text
        postms_dir/
        ├─ Genome/
        |  ├─ Percolator/
        |  └─ post_perc/
        └─ Transcriptome/
           ├─ Percolator/
           └─ post_perc/

    This temp dir mimicks part of the outdir layout for testing, but doesn't
    contain any data.
    """
    base: pathlib.Path = tmp_path_factory.mktemp('postms_dir')
    genome_dir = base / 'Genome'
    transcriptome_dir = base / 'Transcriptome'
    (genome_dir / 'Percolator').mkdir(parents=True)
    (genome_dir / 'post_perc').mkdir(parents=True)
    (transcriptome_dir / 'Percolator').mkdir(parents=True)
    (transcriptome_dir / 'post_perc').mkdir(parents=True)
    return base


@pytest.mark.postms
@pytest.mark.parametrize(
    'folder,filetype',
    [
        ('Genome', 'genome'),
        ('Transcriptome', 'transcriptome')
    ]
)
class TestPostPercolator:
    @pytest.fixture
    def postPerc(
        self,
        data: pathlib.Path,
        folder: str,
        filetype: str,
    ) -> PostPercolator:
        args = Namespace(
            mode='postms',
            genome=data.joinpath("genome.fasta"),
            proteome=data.joinpath("proteome.fasta"),
            gff=data.joinpath("mtb.gtf"),
            mzml=data.joinpath("mzml"),
            rrna=data.joinpath("rrna.fna"),
            transcriptome=True,
        )
        return PostPercolator(args, folder, filetype)

    def test_convert_output(
        self,
        folder: str,
        filetype: str,
        postms_dir: pathlib.Path,
        monkeypatch: pytest.MonkeyPatch,
        data: pathlib.Path,
        postPerc: PostPercolator
    ):
        FILE_NAME = f'{folder}/post_perc/{filetype}_converted_psm.txt'

        # setup
        initial_data = (
            data
            / 'results'
            / folder 
            / 'Percolator'
            / f'{filetype}_results_psm.txt'
        )
        monkeypatch.chdir(postms_dir)
        shutil.copy(initial_data, f'{folder}/Percolator')  # pyright: ignore

        # run
        postPerc.convert_output()
        results = pathlib.Path(FILE_NAME)
        ok_results = data / 'results' / FILE_NAME

        # Assert the expected files exist
        assert results.is_file()
        # Assert resulting data is what we expect
        df = pd.read_csv(results, sep='\t')
        ok_df = pd.read_csv(ok_results, sep='\t')  # pyright: ignore
        pd.testing.assert_frame_equal(df, ok_df, check_like=True)

    # This function tests both the `get_coordinates_genome` and
    # `get_coordinates_rna` from the `PostPercolator` class. Using only one
    # function for both allows us to take advantage of parametrized fixtures,
    # simplifying the code.
    def test_get_coordinates(
        self,
        folder: str,
        filetype: str,
        postms_dir: pathlib.Path,
        monkeypatch: pytest.MonkeyPatch,
        data: pathlib.Path,
        postPerc: PostPercolator
    ):
        FILE_NAME = f'{folder}/post_perc/{filetype}_psm_coords.txt'

        # setup
        initial_data1 = (
            data
            / 'results'
            / folder
            / 'post_perc'
            / f'{filetype}_converted_psm.txt'
        )
        initial_data2 = data / 'assembled.gtf'
        monkeypatch.chdir(postms_dir)
        shutil.copy(initial_data1, f'{folder}/post_perc')
        shutil.copy(initial_data2, '.')

        # run
        if folder == 'Transcriptome': postPerc.get_coordinates_rna()
        else: postPerc.get_coordinates_genome()

        results = pathlib.Path(FILE_NAME)
        ok_results = data / 'results' / FILE_NAME

        # Assert the expected files exist
        assert results.is_file()
        # Assert resulting data is what we expect
        df = pd.read_csv(results, sep='\t')
        ok_df = pd.read_csv(ok_results, sep='\t')  # pyright: ignore
        pd.testing.assert_frame_equal(df, ok_df, check_like=True)

    def test_filter_novel(
        self,
        folder: str,
        filetype: str,
        postms_dir: pathlib.Path,
        monkeypatch: pytest.MonkeyPatch,
        data: pathlib.Path,
        postPerc: PostPercolator
    ):
        FILE_NAME = f'{folder}/post_perc/{filetype}_no_anno.txt'

        # setup
        initial_data = (
            data
            / 'results'
            / folder 
            / 'post_perc'
            / f'{filetype}_psm_coords.txt'
        )
        monkeypatch.chdir(postms_dir)
        shutil.copy(initial_data, f'{folder}/post_perc')  # pyright: ignore

        # run
        postPerc.filter_novel()
        results = pathlib.Path(FILE_NAME)
        ok_results = data / 'results' / FILE_NAME

        # Assert the expected files exist
        assert results.is_file()
        # Assert resulting data is what we expect
        df = pd.read_csv(results, sep='\t')
        ok_df = pd.read_csv(ok_results, sep='\t')  # pyright: ignore
        pd.testing.assert_frame_equal(df, ok_df, check_like=True)

    def test_protein_seqs(
        self,
        folder: str,
        filetype: str,
        postms_dir: pathlib.Path,
        monkeypatch: pytest.MonkeyPatch,
        data: pathlib.Path,
        postPerc: PostPercolator
    ):
        FILE_NAME = f'{folder}/post_perc/{filetype}_proteined.tsv'

        # setup
        initial_data1 = (
            data
            / 'results'
            / folder 
            / 'post_perc'
            / f'{filetype}_linked.tsv'
        )
        initial_data2 = (
            data
            / 'results'
            / f'{filetype}_database.fasta'
        )
        monkeypatch.chdir(postms_dir)
        shutil.copy(initial_data1, f'{folder}/post_perc')  # pyright: ignore
        shutil.copy(initial_data2, '.')  # pyright: ignore

        # run
        postPerc.protein_seqs()
        results = pathlib.Path(FILE_NAME)
        ok_results = data / 'results' / FILE_NAME

        # Assert the expected files exist
        assert results.is_file()
        # Assert resulting data is what we expect
        df = pd.read_csv(results, sep='\t')
        ok_df = pd.read_csv(ok_results, sep='\t')  # pyright: ignore
        pd.testing.assert_frame_equal(df, ok_df, check_like=True)

    def test_add_coordinates(
        self,
        folder: str,
        filetype: str,
        postms_dir: pathlib.Path,
        monkeypatch: pytest.MonkeyPatch,
        data: pathlib.Path,
        postPerc: PostPercolator
    ):
        FILE_NAME = f'{folder}/post_perc/{filetype}_results_02.txt'

        # setup
        initial_data1 = (
            data
            / 'results'
            / folder 
            / 'post_perc'
            / f'{filetype}_proteined.tsv'
        )
        initial_data2 = (
            data
            / 'results'
            / folder 
            / 'post_perc'
            / f'{filetype}_utps.txt'
        )
        monkeypatch.chdir(postms_dir)
        shutil.copy(initial_data1, f'{folder}/post_perc')  # pyright: ignore
        shutil.copy(initial_data2, f'{folder}/post_perc')  # pyright: ignore

        # run
        postPerc.add_coordinates()
        results = pathlib.Path(FILE_NAME)
        ok_results = data / 'results' / FILE_NAME

        # Assert the expected files exist
        assert results.is_file()
        # Assert resulting data is what we expect
        df = pd.read_csv(results, sep='\t')
        ok_df = pd.read_csv(ok_results, sep='\t')  # pyright: ignore
        pd.testing.assert_frame_equal(df, ok_df, check_like=True)
