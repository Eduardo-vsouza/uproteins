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


class Error(Exception):
    """ Base class for other exceptions. """
    pass


class HandlerError(Error):
    """ Raised when the handler is incorrect. """
    def __init__(self, message="Unrecognized handler. Please provide a valid one."):
        self.message = message
        super().__init__(self.message)


class ExternalAssemblyError(Error):
    """ Raised when only one of the two arguments required for using an external transcriptome is provided. """
    def __init__(self, message='Only one of the two arguments required for using an external transcriptome was '
                               'provided. Please, inform both arguments. If using --external_transcriptome, be sure to '
                               'inform --external_gtf as well, and vice-versa.'):
        self.message = message
        super().__init__(self.message)

class FiletypeError(Error):
    """ Raised when the filetype is incorrect. """
    def __init__(self, message='Unrecognized filetype. Please inform one of the following: genome, transcriptome.'):
        self.message = message
        super().__init__(self.message)


class SourceError(Error):
    """ Raised when the source of the data is incorrect. """
    def __init__(self, message="Unrecognized data source. Please provide a valid one."):
        self.message = message
        super().__init__(self.message)


class FormatError(Error):
    """ Raised when the format is incorrect. """
    def __init__(self, message="Unrecognized format. Please provide a valid one."):
        self.message = message
        super().__init__(self.message)


class EngineError(Error):
    """ Raised when the peptide search engine is incorrect. """
    def __init__(self, message="Unsupported peptide search engine. Current supported engines are: percolator, MSGF."):
        self.message = message
        super().__init__(self.message)


class PercolatorProteinsError(Error):
    """ Raised when the file containing proteins from percolator output is missing or incorret."""
    def __init__(self, message="Provide a valid percolator output file containing protein information."):
        self.message = message
        super().__init__(self.message)


class PercolatorPSMError(Error):
    """ Raised when the file containing peptides from percolator output is missing or incorrect."""
    def __init__(self, message="Provide a valid percolator output file containing PSM information."):
        self.message = message
        super().__init__(self.message)

class UProteinsError(Error):
    """ Raised when a file is missing for uProteInS method for identifying unique peptides. """
    def __init__(self, message="Provide both a valid Genbankd and a fasta file containing the predicted ORFs."):
        self.message = message
        super().__init__(self.message)