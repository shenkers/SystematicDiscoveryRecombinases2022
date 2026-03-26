import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import get_junction_reads

get_junction_reads(snakemake.input, snakemake.output)
