import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import kmer_flanks

kmer_flanks(snakemake.input, snakemake.output, snakemake.params)
