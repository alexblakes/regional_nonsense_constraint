import logging
from pathlib import Path
import src
from src.cadd import log_test2

# Module constants
_LOGFILE = f"data/logs/{Path(__file__).stem}.log"

# Logging
logger = logging.getLogger(__name__)

# Functions 
def speak_log(text):
    logger.info(text)

if __name__ == "__main__":
    logger = src.module_logger(_LOGFILE)

    speak_log("Speaking from module 1.")
    log_test2.speak_log("Also speaking from module 1.")
