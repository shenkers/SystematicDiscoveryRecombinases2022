import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import count_raw

count_raw(snakemake.input, snakemake.output, snakemake.params)
