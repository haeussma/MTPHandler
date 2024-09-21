import sys

from loguru import logger

# Flag to check if the logger has already been configured
mtphandler_logger_configured = False


def configure_logger(
    log_level_std: str = "INFO",
    log_level_file: str = "DEBUG",
    log_file: str = "mtp_handler.log",
    log_file_rotation_MB: int = 1,
):
    """Configures the logger severity level for the package."""
    global mtphandler_logger_configured

    if not mtphandler_logger_configured:
        logger.remove()

        logger.add(sys.stdout, level=log_level_std)
        logger.add(
            log_file,
            level=log_level_file,
            rotation=f"{log_file_rotation_MB} MB",
        )
        mtphandler_logger_configured = True
