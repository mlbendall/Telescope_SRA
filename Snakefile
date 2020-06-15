#! /usr/bin/env python
from future import standard_library
standard_library.install_aliases()

import os
import csv
import requests
from subprocess import check_output
from collections import defaultdict


wildcard_constraints:
    proj_acc = "[SED]RP[0-9]{6,9}",
    exp_acc  = "[SED]RX[0-9]{6,9}",
    samp_acc = "[SED]RS[0-9]{6,9}",    
    run_acc  = "[SED]RR[0-9]{6,9}",


configfile: "config.yaml"

"""
Metadata specifications
The metadata table must contain one row for each file. A sample may be comprised as one or
more files are linked by the `sample_id` value.

Important columns include:
  `file_uuid` - File UUID in GDC used to construct download URL
  `md5sum`    - MD5 checksum of the file
  `sample_id` - TCGA format sample id
  `pilot`     - optional, whether this file is part of pilot set
"""

with open(config['sra_run_table'], 'r') as csvfile:
    dialect = csv.Sniffer().sniff(csvfile.read(1024))
    csvfile.seek(0)
    METADATA = [row for row in csv.DictReader(csvfile)]

SAMPLES = [d[config['colname_SAMPLE']] for d in METADATA]

SAMPLE_RUN = defaultdict(list)
for d in METADATA:
    SAMPLE_RUN[d[config['colname_SAMPLE']]].append(d[config['colname_RUN']])


localrules: all
rule all:
    input:
        expand("samples/{s}/completed.txt", s=SAMPLES)


include: "refs.smk"
include: "runs.smk"
include: "samples.smk"
