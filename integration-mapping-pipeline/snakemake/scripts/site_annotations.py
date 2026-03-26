import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import get_site_annotations

get_site_annotations(snakemake.input, snakemake.output, snakemake.params)
