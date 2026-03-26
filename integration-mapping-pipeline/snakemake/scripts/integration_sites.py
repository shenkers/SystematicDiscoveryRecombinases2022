import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import get_integration_sites

get_integration_sites(snakemake.input, snakemake.output, snakemake.params)
