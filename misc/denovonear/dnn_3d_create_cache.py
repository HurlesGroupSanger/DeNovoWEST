#!/usr/bin/env python

import asyncio

import click
from denovonear.load_de_novos import load_de_novos
from denovonear.load_gene import (
    get_uniprot_ids_for_transcript,
    load_gene,
)
from denovonear.rate_limiter import RateLimiter
from gencodegenes.gencode import Gencode


async def _load_gencode(symbols, build="grch37"):
    """load gene coords and sequence via ensembl"""
    gencode = Gencode()
    async with RateLimiter(per_second=15) as ensembl:
        tasks = [load_gene(ensembl, symbol, build) for symbol in symbols]
        genes = await asyncio.gather(*tasks)
        for gene in genes:
            gencode.add_gene(gene)
        return gencode


def load_gencode(symbols, gencode=None, fasta=None, build="grch37"):
    """load genes from gencode annotations file, with ensembl as backup"""
    if gencode and fasta:
        return Gencode(gencode, fasta)

    # use ensembl as backup if gencode file not available. This restricts the
    # asynchronous calls to within one section, and ensures they are called together
    return asyncio.get_event_loop().run_until_complete(_load_gencode(symbols, build))


@click.command()
@click.argument("dnm")
@click.argument("gencode_gtf")
@click.argument("fasta")
@click.argument("build", default="grch38", type=click.Choice(["grch37", "grch38"]))
def main(dnm, gencode_gtf, fasta, build):
    """
    When running denovonear-3D on a HPC, there are issues of concurrent writing access to the Ensembl cache.
    We workaround this issue by running a prelimimary process that builds the cache.

    Args:
        dnm (str): path to de novo mutations file
        gencode_gtf (str): path to gencode gtf file
        fasta (str): path to fasta file
        build (str): genome build, either grch37 or grch38
    """

    de_novos = load_de_novos(dnm)
    gencode = load_gencode(de_novos, gencode_gtf, fasta, build)

    cpt = 0
    per_gene_uniprot_ids = {}
    for symbol in sorted(de_novos):
        if len(de_novos[symbol]["missense"] + de_novos[symbol]["nonsense"]) < 2:
            continue

        if symbol not in gencode:
            continue

        uniprot_ids = get_uniprot_ids_for_transcript(gencode[symbol].canonical.name, "grch38")
        per_gene_uniprot_ids[symbol] = uniprot_ids
        cpt += 1
        if cpt % 100 == 0:
            print(cpt)

    with open("ensembl_cache_gene_uniprot_ids.txt", "w") as f:
        for symbol, uniprot_ids in per_gene_uniprot_ids.items():
            f.write(f"{symbol}:{','.join(uniprot_ids)}\n")


if __name__ == "__main__":
    main()
