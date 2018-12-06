# -*- coding: utf-8 -*-

"""
manorm.logging
~~~~~~~~~~~~~~

Logging configurations for MAnorm.
"""

from __future__ import absolute_import

import logging
import sys

# Initialize the default logger with a NullHandler
logger = logging.getLogger('manorm')
logger.addHandler(logging.NullHandler())


class CleanFormatter(logging.Formatter):
    """Clean logging formatter, omit the level label for logging.INFO."""

    def format(self, record):
        if record.levelno != logging.INFO:
            record.msg = "{}: {}".format(record.levelname, record.msg)
        return super(CleanFormatter, self).format(record)


def setup_logger(verbose=False):
    """Setup MAnorm logger with a stderr stream handler under given verbose level."""
    logger.setLevel(logging.DEBUG)
    # clear all handlers
    for handler in logger.handlers:
        logger.removeHandler(handler)

    sh = logging.StreamHandler(stream=sys.stderr)
    if verbose:
        sh.setLevel(logging.DEBUG)
        formatter = logging.Formatter("%(asctime)s %(name)-20s %(lineno)-4d %(levelname)-8s %(message)s",
                                      datefmt="%Y-%m-%d %H:%M")
    else:
        sh.setLevel(logging.INFO)
        formatter = CleanFormatter()
    sh.setFormatter(formatter)
    logger.addHandler(sh)
