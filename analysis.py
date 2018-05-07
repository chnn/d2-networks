# See README for instructions.

import json
from glob import glob
from itertools import product
from os import getenv

from d2 import d2
from phylum_data import PHYLUM_DATA
from pyfasta import Fasta

K = 25

seq_data = {}
scores = {}
metadata = {}

i = 0

for filename in glob(getenv("DATA_DIR", "data") + "/*.fna"):
    fasta = Fasta(filename)
    key = sorted(fasta.keys())[0]

    genbank_id = key.split(" ")[0]
    short_name = " ".join(key.split(" ")[1:3])
    org_phylum_data = PHYLUM_DATA.get(short_name, {})
    name = " ".join(key.split(" ")[1:-2])[:-1]

    metadata[genbank_id] = {
        "name": name,
        "phylum": org_phylum_data.get("phylum", ""),
        "domain": org_phylum_data.get("domain", ""),
        "ncbiLevel3": org_phylum_data.get("ncbiLevel3", "")
    }

    seq_data[genbank_id] = fasta[key][:]

for k1, k2 in product(seq_data.keys(), seq_data.keys()):
    if k1 == k2:
        continue

    key = "({}, {})".format(*sorted([k1, k2]))

    if key not in scores:
        scores[key] = d2(seq_data[k1], seq_data[k2], K)

print(json.dumps({"scores": scores, "metadata": metadata}))
