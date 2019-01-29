#!/usr/bin/env python3
#
# MIDAS IGG SMELTER
#
# A tool to repackage and enrich information provided by Stephen Nayfach
# via http://github.com/snayfach/IGGdb in order to support genotyping and
# SNP calling over the IGGdb dataset.
#
# Created on January 28, 2019 by Boris Dimitrov, bdimitrov@chanzuckerberg.com.
# Distributed freely under the MIT License, subject to restrictions set forth
# in the licensing terms of MIDAS and IGG tools by Stephen Nayfach.  Please
# see http://github.com/snayfach/MIDAS and http://github.com/snayfach/IGGsearch
# for appropriate citation requirements.

import sys
assert sys.version_info >= (3, 6), "Please run this script with Python version >= 3.6."


from os.path import dirname, basename, abspath, isdir
import traceback
import json
import random
import time
from utilities import tsprint, backtick, makedirs, tsv_rows, parse_table, ProgressTracker


def smelt(argv):
    cwd = backtick('pwd')
    tsprint(
        f"cd {cwd}\n" +
        ' '.join(argv)
    )
    _, subcmd, outdir, iggtoc = argv
    assert subcmd.replace("-", "_").lower() == "build_gsnap_index"
    assert basename(iggtoc) == "species_info.tsv"
    try:
        metadata = dirname(abspath(iggtoc))
        assert basename(metadata) == "metadata"
        iggdb = dirname(metadata)
        # assert re.match("^v[0-9]+\.[0-9]+\.[0-9]+$", basename(iggdb))
        pangenomes = f"{iggdb}/pangenomes"
        assert isdir(pangenomes)
    except Exception as e:
        e.help_text = f"Unexpected directory structure for MIDAS-IGGdb database around {iggtoc}."
        raise
    makedirs(outdir, exist_ok=False)
    species = list(parse_table(tsv_rows(iggtoc)))
    tsprint(f"Found {len(species)} species in {iggtoc}, for example:")
    random.seed(time.time())
    random_index = random.randrange(0, len(species))
    tsprint(json.dumps(species[random_index], indent=4))
    tsprint("Now collating fasta for gsnap index construction.")
    MAX_FAILURES = 100
    count_successes = 0
    ticker = ProgressTracker(target=len(species))
    failures = []
    for s in species:
        # Prepend species_alt_id to each header in the species pangenome, and append to all.fa.
        try:
            s_species_alt_id = s['species_alt_id']
            s_pangenome = f"{pangenomes}/{s_species_alt_id}/centroids.fa"
            s_tempfile = f"{outdir}/temp_{s_species_alt_id}.fa"
            s_header_xform = f"sed 's=^>=>{s_species_alt_id}|=' {s_pangenome} > {s_tempfile} && cat {s_tempfile} >> {outdir}/temp_pangenomes.fa && rm {s_tempfile} && echo SUCCEEDED || echo FAILED"
            status = backtick(s_header_xform)
            assert status == "SUCCEEDED"
            count_successes += 1
        except Exception as e:
            failures.append(s)
            if len(failures) == MAX_FAILURES:
                count_examined = len(failures) + count_successes
                e.help_text = f"Giving up after {MAX_FAILURES} failures in first {count_examined} species.  See temp files for more info."
                raise
        finally:
            ticker.advance(1)
    failed_species_alt_ids = [s['species_alt_id'] for s in failures]
    if not failures:
        tsprint(f"All {len(species)} species were processed successfully.")
    else:
        tsprint(f"Collation of {len(failures)} species failed.  Those are missing from the final pangenomes.fa")
    # Create output file only on success.
    # Dump stats in json.
    collation_status = {
        "comment": f"Collation into pangenomes.fa succeeded on {time.asctime()}.",
        "successfully_collated_species_count": count_successes,
        "failed_species_count": len(failures),
        "total_species_count": len(species),
        "failed_species_alt_ids": failed_species_alt_ids
    }
    collation_status_str = json.dumps(collation_status, indent=4)
    with open(f"{outdir}/pangenomes_collation_status.json", "w") as pcs:
        chars_written = pcs.write(collation_status_str)
        assert chars_written == len(collation_status_str)
        tsprint(collation_status_str)
    rename(f"{outdir}/temp_pangenomes.fa", f"{outdir}/pangenomes.fa")

def main():
    try:
        smelt(sys.argv)
    except Exception as e:
        tsprint(traceback.format_exc())
        tsprint("*** USAGE:  See https://github.com/czbiohub/MIDAS-IGGdb/blob/master/README.md#smelter ***\n")
        if hasattr(e, 'help_text'):
            tsprint(f"*** {e.help_text} ***") # pylint: disable=no-member


if __name__ == "__main__":
    main()
