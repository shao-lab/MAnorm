import logging

from manorm.logging import logger, setup_logger, CleanFormatter


def test_clean_formatter():
    formatter = CleanFormatter()
    record = logging.makeLogRecord(
        dict(levelname="INFO", levelno=logging.INFO, msg='message info'))
    assert formatter.format(record) == "message info"
    record = logging.makeLogRecord(
        dict(levelname="DEBUG", levelno=logging.DEBUG, msg='message debug'))
    assert formatter.format(record) == "DEBUG: message debug"


def test_setup_logger():
    assert len(logger.handlers) == 1
    assert isinstance(logger.handlers[0], logging.NullHandler)
    setup_logger(verbose=False)
    assert logger.level == logging.DEBUG
    assert len(logger.handlers) == 1
    assert isinstance(logger.handlers[0], logging.StreamHandler)
    assert logger.handlers[0].level == logging.INFO
    setup_logger(verbose=True)
    assert logger.level == logging.DEBUG
    assert len(logger.handlers) == 1
    assert isinstance(logger.handlers[0], logging.StreamHandler)
    assert logger.handlers[0].level == logging.DEBUG
