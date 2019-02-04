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


import os
from os.path import dirname, basename, abspath, isdir
import traceback
import json
import random
import time
from utilities import tsprint, backtick, makedirs, ProgressTracker
from iggdb import IGGdb


def smelt(argv):
    cwd = backtick('pwd')
    tsprint(
        f"cd {cwd}\n" +
        ' '.join(argv)
    )
    _, subcmd, outdir, iggtoc = argv
    subcmd = subcmd.replace("-", "_").lower()
    SUBCOMMANDS = {
        f"build_gsnap_{gdim}_index": gdim
        for gdim in ["pangenomes", "repgenomes"]
    }
    gdim = SUBCOMMANDS.get(subcmd)
    try:
        assert gdim in ["pangenomes", "repgenomes"]
    except Exception as e:
        e.help_text = f"Try a supported subcommand instead of {subcmd}."
        raise
    assert basename(iggtoc) == "species_info.tsv"
    try:
        metadata = dirname(abspath(iggtoc))
        assert basename(metadata) == "metadata"
        iggdb_root = dirname(metadata)
        assert isdir(f"{iggdb_root}/pangenomes")
        assert isdir(f"{iggdb_root}/repgenomes")
    except Exception as e:
        e.help_text = f"Unexpected directory structure for MIDAS-IGGdb database around {iggtoc}."
        raise
    makedirs(outdir, exist_ok=False)
    iggdb = IGGdb(iggdb_root)
    tsprint(f"Found {len(iggdb.species)} species in {iggtoc}, for example:")
    random.seed(time.time())
    random_index = random.randrange(0, len(iggdb.species))
    tsprint(json.dumps(iggdb.species[random_index], indent=4))
    tsprint(f"Now collating fasta for gsnap {gdim} index construction.")
    MAX_FAILURES = 100
    count_successes = 0
    ticker = ProgressTracker(target=len(iggdb.species))
    failures = []
    for s in iggdb.species:
        try:
            s_species_alt_id = s['species_alt_id']
            s_tempfile = f"{outdir}/temp_{gdim}_{s_species_alt_id}.fa"
            # Note how the header tags we emit below for pangenomes and repgenomes are consistent;
            # this should enable easy reconciliation of gsnap alignments against the
            # two separate indexes.
            if gdim == "pangenomes":
                #
                # The header tag we wish to emit would be
                #
                #    >4547837|1657.8.patric|rest_of_original_header_from_pangenome_file
                #
                # where
                #
                #    species_alt_id = 4547837                                # from species_alt_id column in table
                #    s_rg_id_origin = 1657.8.patric                          # from original header in pangenome file
                #
                # As the original header in the pangenome file already begins
                # with s_rg_id_origin, we just need to prepend species_alt_id.
                #
                s_pangenome = f"{iggdb_root}/pangenomes/{s_species_alt_id}/centroids.fa"
                s_header_xform = f"sed 's=^>=>{s_species_alt_id}|=' {s_pangenome} > {s_tempfile} && cat {s_tempfile} >> {outdir}/temp_{gdim}.fa && rm {s_tempfile} && echo SUCCEEDED || echo FAILED"
            else:
                assert gdim == "repgenomes"
                #
                # The header tag we wish to emit would be
                #
                #    >4547837|1657.8.patric|entire_original_header_from_repgenome_file
                #
                # where
                #
                #    species_alt_id = 4547837                                # species_alt_id column in species_info
                #    s_rg_id_origin = 1657.8.patric                          # from file listing in repgenomes dir
                #
                s_rg_id_origin = s['representative_genome_with_origin']
                s_rg_path = s['representative_genome_path']
                s_header_xform = f"sed 's=^>=>{s_species_alt_id}|{s_rg_id_origin}|=' {s_rg_path} > {s_tempfile} && cat {s_tempfile} >> {outdir}/temp_{gdim}.fa && rm {s_tempfile} && echo SUCCEEDED || echo FAILED"
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
        tsprint(f"All {len(iggdb.species)} species were processed successfully.")
    else:
        tsprint(f"Collation of {len(failures)} species failed.  Those are missing from the final {gdim}.fa")
    # Create output file only on success.
    # Dump stats in json.
    collation_status = {
        "comment": f"Collation into {gdim}.fa succeeded on {time.asctime()}.",
        "successfully_collated_species_count": count_successes,
        "failed_species_count": len(failures),
        "total_species_count": len(iggdb.species),
        "failed_species_alt_ids": failed_species_alt_ids
    }
    collation_status_str = json.dumps(collation_status, indent=4)
    with open(f"{outdir}/{gdim}_collation_status.json", "w") as pcs:
        chars_written = pcs.write(collation_status_str)
        assert chars_written == len(collation_status_str)
        tsprint(collation_status_str)
    os.rename(f"{outdir}/temp_{gdim}.fa", f"{outdir}/{gdim}.fa")

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