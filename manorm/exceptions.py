"""
manorm.exceptions
-----------------

This module contains the exceptions of MAnorm.
"""


class ManormError(Exception):
    """Base class for MAnorm errors."""


class FileFormatError(ManormError):
    def __init__(self, format, line_num, line):
        msg = f"Invalid {format} format at line {line_num}: {line!r}"
        super().__init__(msg)


class FormatModeConflictError(Exception):
    def __init__(self, format, mode):
        msg = f"format {format!r} is conflict with mode {mode!r}"
        super().__init__(msg)


class ProcessNotReadyError(ManormError):
    def __init__(self, step, precursor):
        msg = f"Unable to {step}, please {precursor} first"
        super().__init__(msg)
