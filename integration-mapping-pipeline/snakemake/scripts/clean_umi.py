import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import clean_umi

clean_umi(snakemake.input, snakemake.output, snakemake.params)
