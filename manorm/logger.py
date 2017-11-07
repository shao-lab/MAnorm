""""Configures the logger of MAnorm."""

import sys
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
ch = logging.StreamHandler(stream=sys.stderr)
ch.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)-8s %(asctime)s: %(message)s", datefmt="%m/%d/%Y %H:%M")
ch.setFormatter(formatter)
logger.addHandler(ch)
