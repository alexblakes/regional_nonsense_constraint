import logging
from pathlib import Path
import src
from src.cadd import log_test1

# Module constants
_LOGFILE = f"data/logs/{Path(__file__).stem}.log"

# Logging
logger = logging.getLogger(__name__)

# Functions 
def speak_log(text):
    logger.info(text)

if __name__ == "__main__":
    logger = logging.getLogger("src")
    src.add_filehandler(logger, _LOGFILE, mode="w")

    speak_log("Speaking from module 2.")
    log_test1.speak_log("Also speaking from module 2.")