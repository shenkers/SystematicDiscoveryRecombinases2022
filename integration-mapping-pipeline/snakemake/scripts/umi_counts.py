import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import count_umis

count_umis(snakemake.input, snakemake.output, snakemake.params)
