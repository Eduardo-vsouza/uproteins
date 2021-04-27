class Error(Exception):
    """ Base class for other exceptions. """
    pass


class HandlerError(Error):
    """ Raised when the handler is incorrect. """
    def __init__(self, message="Unrecognized handler. Please provide a valid one."):
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