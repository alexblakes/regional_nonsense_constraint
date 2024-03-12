""""Package docstring."""

import logging

from src import constants as C

# Functions
def log_to_file():
    file_handler = logging.FileHandler(f"data/logs/testing.log", "w")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)        

    logging.addHandler(file_handler)


# Logging configuration
formatter = logging.Formatter(
    fmt="[%(asctime)s] %(levelname)s %(funcName)s(): %(message)s", 
    datefmt='%d-%b-%y %H:%M:%S',
    )

file_handler = logging.FileHandler(f"data/logs/testing.log", "w")
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(formatter)        

stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.INFO)
stream_handler.setFormatter(formatter)

root = logging.getLogger()
root.setLevel(logging.DEBUG)
root.addHandler(stream_handler)
root.addHandler(file_handler)