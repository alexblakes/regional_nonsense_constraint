import logging
from pathlib import Path
import src

# Module constants
_LOGFILE = f"data/logs/{Path(__file__).stem}.log"

# Logging
logger = logging.getLogger(__name__)

# Functions

# __main__ clause
if __name__ == "__main__":
    logger = src.module_logger(_LOGFILE)
