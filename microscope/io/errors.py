"""Shared I/O exceptions."""


class FileFormatError(Exception):
    """The file could not be parsed as the expected format."""


class UnsupportedFormatError(Exception):
    """No parser is available for this file."""
