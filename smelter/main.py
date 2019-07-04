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
import traceback
import json
import time
import threading
import multiprocessing
from collections import defaultdict
from utilities import tsprint, backtick, run, makedirs, ProgressTracker, tsv_rows, parse_table

from iggdb import IGGdb


def smelt(argv):
    cwd = backtick('pwd')
    my_command = f"cd {cwd}; " + ' '.join(argv)
    tsprint(my_command)
    _, subcmd, outdir, iggdb_toc = argv
    subcmd = subcmd.replace("-", "_").lower()
    SUBCOMMANDS = {
        f"collate_{gdim}": (gdim, collate)
        for gdim in ["pangenomes", "repgenomes"]
    }
    SUBCOMMANDS.update({
        "rename_markers": (None, rename_markers)
    })
    try:
        gdim, gfunc = SUBCOMMANDS[subcmd]
    except Exception as e:
        e.help_text = f"Try a supported subcommand instead of {subcmd}.  Supported subcommands: {SUBCOMMANDS.keys()}."
        raise
    makedirs(outdir, exist_ok=False)
    iggdb = IGGdb(iggdb_toc)
    gfunc(iggdb, gdim, outdir, my_command)


def collate(iggdb, gdim, outdir, my_command):
    tsprint(f"Now collating fasta for gsnap {gdim} index construction.")
    MAX_FAILURES = 100
    count_successes = 0
    ticker = ProgressTracker(target=len(iggdb.species_info))
    failures = []
    for s in iggdb.species_info:
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
                #    species_alt_id = 4547837                  # from species_alt_id column in table
                #    repgenome_with_origin = 1657.8.patric     # from original header in pangenome file
                #
                # As the original header in the pangenome file already begins
                # with s_rg_id_origin, we just need to prepend species_alt_id.
                #
                s_header_xform = f"sed 's=^>=>{s_species_alt_id}|=' {s['pangenome_path']} > {s_tempfile} && cat {s_tempfile} >> {outdir}/temp_{gdim}.fa && rm {s_tempfile} && echo SUCCEEDED || echo FAILED"
            else:
                assert gdim == "repgenomes"
                #
                # The header tag we wish to emit would be
                #
                #    >4547837|1657.8.patric|entire_original_header_from_repgenome_file
                #
                # where
                #
                #    species_alt_id = 4547837                # species_alt_id column in species_info
                #    repgenome_with_origin = 1657.8.patric   # from file listing in repgenomes dir
                #
                s_repgenome_with_origin = s['repgenome_with_origin']
                s_repgenome_path = s['repgenome_path']
                s_header_xform = f"sed 's=^>=>{s_species_alt_id}|{s_repgenome_with_origin}|=' {s_repgenome_path} > {s_tempfile} && cat {s_tempfile} >> {outdir}/temp_{gdim}.fa && rm {s_tempfile} && echo SUCCEEDED || echo FAILED"
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
        tsprint(f"All {len(iggdb.species_info)} species were processed successfully.")
    else:
        tsprint(f"Collation of {len(failures)} species failed.  Those are missing from the final {gdim}.fa")
    # Create output file only on success.
    # Dump stats in json.
    collation_status = {
        "comment": f"Collation into {gdim}.fa succeeded on {time.asctime()} with command '{my_command}'.",
        "successfully_collated_species_count": count_successes,
        "failed_species_count": len(failures),
        "total_species_count": len(iggdb.species_info),
        "failed_species_alt_ids": failed_species_alt_ids,
        "elapsed_time": time.time() - ticker.t_start
    }
    collation_status_str = json.dumps(collation_status, indent=4)
    with open(f"{outdir}/{gdim}_collation_status.json", "w") as pcs:
        chars_written = pcs.write(collation_status_str)
        assert chars_written == len(collation_status_str)
        tsprint(collation_status_str)
    os.rename(f"{outdir}/temp_{gdim}.fa", f"{outdir}/{gdim}.fa")


def fetch_remap(spid, precgenes, MAX_ATTEMPTS=2):
    for attempt in range(1, MAX_ATTEMPTS + 1):
        result = []
        try:
            blast_output = backtick(f"aws s3 cp  s3://microbiome-chunyu/midas-iggdb/dbs/prokka_to_pan/{spid}.blastn.tsv -")
        except:
            if attempt == MAX_ATTEMPTS:
                raise
            tsprint(f"WARNING:  Retrying aws cp for secies {spid}.")
            time.sleep(10)
    for lno, line in enumerate(blast_output.split("\n")):
        try:
            prokka_id, iggdb_id = line.split()[:2]
            if prokka_id in precgenes:
                # tsprint(f"{prokka_id} -> {iggdb_id}")
                result.append((prokka_id, iggdb_id))
        except:
            tsprint(f"ERROR:{spid}:{lno}: " + line)
            raise
    return result


def rename_markers(iggdb, gdim, outdir, _my_command):
    assert gdim == None
    iggdb_root = iggdb.iggdb_root
    # from s3://microbiome-chunyu/midas-iggdb/dbs/marker_genes
    # marker genes generated by chunyu by running hmmsearch on the 23,790 repgenomes of iggdb
    # (even though we need to run on all 206,000+ genomes, we are starting with the 23k repgenomes)
    # "precursor" not in a biological sense but because gene names are autogenerated by prokka and
    # do not match the naming scheme for IGGdb/midas;  to be translated from prokka to iggdb
    # naming scheme by the artisanally handcrafted code below
    precursor_file = f"{iggdb_root}/metadata/precursor_marker_genes/phyeco.map"
    precursors = list(parse_table(tsv_rows(precursor_file)))
    precursor_genes = defaultdict(set)
    precursor_markers = defaultdict(set)
    species = set(s['species_alt_id'] for s in iggdb.species.values())
    hist = defaultdict(int)
    for p in precursors:
        assert p['species_id'] in species, f"Unlisted species {p['species_id']}.  Please add to {iggdb.iggdb_toc_species}."
        precursor_genes[p['species_id']].add(p['gene_id'])
        precursor_markers[p['species_id']].add(p['marker_id'])
    counts_str = "\t".join(("species_id", "count_marker_genes", "count_marker_gene_families", "multimapped(!)")) + "\n"
    for s in sorted(species):
        multimapped = "MULTIMAPPED!" if len(precursor_genes[s]) != len(precursor_markers[s]) else ""
        counts_str += "\t".join(str(val) for val in (s, len(precursor_genes[s]), len(precursor_markers[s]), multimapped)) + "\n"
        hist[len(precursor_markers[s])] += 1
    summary_file = f"{outdir}/precursor_phyeco_summary.tsv"
    with open(summary_file, "w") as pps:
        pps.write(counts_str)
    hist_str = "\t".join(("marker_gene_families_count", "species_count", "cumulative_percent")) + "\n"
    sofar = 0
    total = len(species)
    for marker_gene_families_count, species_count in sorted(hist.items()):
        sofar += species_count
        hist_str += f"{marker_gene_families_count}\t{species_count}\t{100.0*sofar/total:3.1f}%\n"
    hist_file = f"{outdir}/precursor_phyeco_histogram.tsv"
    with open(hist_file, "w") as hf:
        hf.write(hist_str)
    run(f"cat {hist_file}")
    mp = multiprocessing.Pool(2 * multiprocessing.cpu_count())
    remap = {}
    t_start = time.time()
    work = sorted(species)
    for submap in mp.starmap(fetch_remap, ((spid, precursor_genes[spid]) for spid in work)):
        for k, v in submap:
            remap[k] = v
    t_end = time.time()
    tsprint(f"Fetched {len(work)} remaps in {t_end - t_start:3.1f} seconds.")
    missing_remaps_genes = set()
    missing_remaps_species = set()
    for p in precursors:
        p["prokka_id"] = p["gene_id"]
        p["gene_id"] = remap.get(p["prokka_id"])
        if not p["gene_id"]:
            missing_remaps_genes.add(p["prokka_id"])
            missing_remaps_species.add(p["species_id"])
        p["species_alt_id"] = p["species_id"]
        del p["species_id"]
        del p["genome_id"]
    if missing_remaps_genes:
        tsprint(f"Missing remaps for {len(missing_remaps_genes)} genes from {len(missing_remaps_species)} species.")
        assert len(missing_remaps_genes) <= 0.1 * len(precursors)
        assert len(missing_remaps_species) <= 0.1 * len(species)
    cols = list(precursors[0].keys())
    remapped_str = "\t".join(cols) + "\n"
    tsprint(remapped_str)
    for p in precursors:
        try:
            remapped_str += "\t".join(p[c] or "None" for c in cols) + "\n"
        except:
            tsprint(json.dumps(p, indent=4))
            raise
    remapped_file = f"{outdir}/phyeco.map"
    with open(remapped_file, "w") as rf:
        rf.write(remapped_str)
    run(f"head {remapped_file}")


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
