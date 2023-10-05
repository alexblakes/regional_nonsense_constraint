""""Module docstring."""

import logging

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s [%(module)s;%(funcName)s] %(message)s",
    datefmt='%d-%b-%y %H:%M:%S'
)
