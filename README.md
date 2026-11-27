# CarveMe + SMETANA on TheBelly

This package takes one **protein FASTA per major bin**, reconstructs metabolic
models with CarveMe, and runs SMETANA on the recovered members of each LT
culture. It uses CPU resources. The defaults are four concurrent reconstructions
with eight DIAMOND threads each; SMETANA analyses run sequentially.

The default solver is **SCIP**, through PySCIPOpt. No commercial solver licence
is needed. The supplied installer requires an existing conda or mamba installation
and creates a new environment without sudo.

## 1. Install once

```bash
unzip lt_carveme_smetana.zip
cd lt_carveme_smetana
bash install.sh
conda activate lt-metabolism
```

If conda is not installed, install [Miniforge](https://github.com/conda-forge/miniforge)
for Linux x86_64 first. If `conda activate` is unavailable in a shell, initialize
the conda shell hook according to that installation's instructions.

The installer pins CarveMe 1.6.6, SMETANA 1.2.1, REFRAMED 1.6.0,
PySCIPOpt 6.3.0, and DIAMOND 2.1.11, using Python 3.12. It then runs a real
SBML/FBA/SMETANA cross-feeding test. A successful installation ends with `PASS`.
The CarveMe wheel includes its standard bacterial universe and reference data.
The resolved Python package versions and DIAMOND version are recorded with
each production run.

## 2. Prepare the inputs

Use the per-bin `.faa` files from gene prediction or your annotation workflow.
Existing eggNOG annotations can be used separately to assess the predicted
pathways; this pipeline starts from the protein sequences and does not require
eggNOG-mapper to be run again. Do not supply nucleotide contigs as `.faa`.

The automatic mode reads `.faa` or `.faa.gz` directly inside one directory.
Filenames must contain the sample and genus, for example:

```text
LT1-Geobacter.faa
LT2-Geobacter.faa
LT3-Geobacter.faa
LT3-Azospira.faa
LT4-Geobacter.faa
LT4-Azospira.faa
LT4-Cupidesulfovibrio.faa
LT5-Geobacter.faa
LT5-Azospira.faa
LT5-Cupidesulfovibrio.faa
```

Your previous `Nitratidesulfovibrio` filenames are accepted and mapped to
`Cupidesulfovibrio` in the output, following the revised phylogeny. Original
files are preserved. Files for the small Oryzomicrobium, Cellulomonas, and
Rectinema bins are excluded from automatic selection. Their incompleteness
would make missing pathways especially difficult to interpret.

If the filenames differ, copy `manifest.example.tsv`, replace its paths, and
use `--manifest my_bins.tsv` in place of `--proteins /path/to/proteins`.
Paths in the manifest are relative to the manifest, or can be absolute.
The ten expected major bins are required; duplicate sample/genus matches stop
the run so that a second assembly is not silently included.

Check the selection before starting any modelling:

```bash
python lt_metabolism.py \
    --proteins /absolute/path/to/bin_proteins \
    --outdir /absolute/path/to/LT_metabolism \
    --plan-only
```

If you only have contig FASTAs, predict proteins first. For example, with
Prodigal installed, `prodigal -i LT4-Geobacter.fna -a LT4-Geobacter.faa -p single`.
Do this separately for each bin, preserving the sample/genus filenames.

## 3. First pass: reconstruct and explore

```bash
python lt_metabolism.py \
    --proteins /absolute/path/to/bin_proteins \
    --outdir /absolute/path/to/LT_metabolism \
    --jobs 4 --threads 8 \
    --exploratory --pairs
```

This builds ten models. The observed multispecies groups are:

| Culture | Models included |
|---|---|
| LT1, LT2 | Geobacter only; reconstructed and checked for growth |
| LT3 | Geobacter + Azospira |
| LT4 | Geobacter + Azospira + Cupidesulfovibrio |
| LT5 | Geobacter + Azospira + Cupidesulfovibrio |

`--pairs` also tests the three two-member subsets of LT4 and LT5. They are
model comparisons, rather than additional observed cultures. No analysis
combines organisms from different LT samples.

**Exploratory conditions:** global MIP/MRO uses SMETANA's complete anaerobic
environment; detailed SMETANA uses its inferred anaerobic minimal environment.
The medium labels are retained in the output. These two exploratory modes
therefore describe different environments and should not be directly compared
as if they used the same recipe. Their purpose is to check the reconstruction
and identify candidates before defining the experimental medium.

To reconstruct models without running SMETANA, substitute `--models-only`
for `--exploratory --pairs`.

## 4. Run the actual assay medium

For the paper, encode the medium used in the original physiological assay.
The culture maintenance medium and the phenotyping medium may differ; run
separate conditions where appropriate.

SMETANA's medium database is a two-column, tab-separated file with a header:

```text
medium	compound
assay	ac
assay	fum
assay	nh4
assay	pi
assay	h2o
```

**Those five rows illustrate the format; they are not a complete recipe.**
Make `assay_media.tsv` from the actual recipe, including inorganic nutrients,
trace nutrients, vitamins, donor, and acceptor as represented by the models.
Use bare BiGG metabolite identifiers: `ac`, `fum`, `for`, `h2`, `so4`, etc.,
without `M_`, `EX_`, or compartment suffixes. Oxygen is rejected in a selected
medium because this pipeline is for anaerobic culture conditions.

```bash
python lt_metabolism.py \
    --proteins /absolute/path/to/bin_proteins \
    --outdir /absolute/path/to/LT_metabolism \
    --resume --pairs \
    --media-db /absolute/path/to/assay_media.tsv \
    --media assay
```

Multiple defined recipes can be compared with `--media maintenance,assay`.
Each selected recipe is applied to all cultures for a consistent comparison.
These medium files specify availability, **not measured concentrations or
uptake rates**; SMETANA uses its default uptake limits. Neither `--anaerobic`
alone nor a filename establishes the experimental medium.

`--resume` reuses models only if input hashes, software versions, settings,
and model checksums match. Each analysis gets a new dated directory, so a
second medium can be analysed without replacing earlier results. Use a new
output directory for changed proteins or a changed reconstruction environment.
Do not run two copies against the same output directory simultaneously.

## 5. Read the results

```text
LT_metabolism/
  resolved_manifest.tsv
  inputs/                         normalized proteins and original-header maps
  models/                         ten SBML models and their provenance JSON
  logs/                           CarveMe logs and installed versions
  work/                           reconstruction intermediates
  analysis/<dated_run>/
    settings.json
    media.tsv                     present for a supplied medium
    monoculture_growth.tsv
    model_exchanges.tsv
    medium_compounds_without_exchange.tsv
    run_status.tsv
    all_global.tsv
    all_detailed.tsv
    positive_exchanges.tsv
    LT3/, LT4/, LT5/, ...          per-community results, commands, and logs
```

Start with `run_status.tsv` and `monoculture_growth.tsv`. The pipeline checks
whether each community can support every member at SMETANA's default biomass
threshold of 0.1 in that mode's environment. An infeasible community is recorded
as `skipped_no_joint_growth`, rather than assigned a score of zero. Individual
models can fail to grow while the corresponding community grows; that can be
a useful candidate dependency, provided the missing pathway is credible.
ATP maintenance is set to zero in the diagnostic growth calculations to match
SMETANA's own model-loading behaviour.

`positive_exchanges.tsv` contains predicted donor, receiver, metabolite, and
SMETANA score. It is sorted by score. The score combines species coupling
(SCS), uptake requirement (MUS), and production potential (MPS). It is not a
flux estimate or an experimentally calibrated probability of interaction.
All detailed rows, including zeros, are retained in `all_detailed.tsv`.
`review_log` marks solver/algorithm warnings or missing global scores even
when the underlying command exited successfully. A missing result and a
computed zero require different interpretations.

Global outputs report metabolic interaction potential (MIP) and resource
overlap (MRO). Use the detailed results to identify specific exchanges rather
than inferring metabolite transfer from a global score alone.

For your study, inspect **both directions** between Geobacter and each other
organism. Hydrogen, formate, organic acids, sulfur compounds, and growth
factors are reasonable candidates to examine when present in the reconstructed
networks. Compare the pairs with the three-member cultures, then inspect the
supporting reactions, transporters, and genes in the relevant bins. The gene-ID
maps link reconstructed gene identifiers back to your original protein headers.

CarveMe's standard bacterial universe is a starting point. Extracellular
electron transfer to an electrode or insoluble Fe(III), direct interspecies
electron transfer, and abiotic reactions of sulfide with Fe(III) require
additional curation or experiments. Adding an iron salt to a medium list does
not establish a correctly modelled respiratory pathway. Availability of the
terminal electron acceptor is especially important before explaining a
Geobacter growth phenotype with these models.

No additional medium-specific gap filling is requested by this pipeline.
CarveMe still performs its normal optimisation and can retain unsupported
reactions to build a viable model. Filling every organism to grow independently
on the assay medium could erase the dependencies being investigated, so review
any proposed manual gap filling alongside gene evidence.

## A compatibility correction included here

During testing, unmodified SMETANA 1.2.1 + REFRAMED 1.6.0 missed an obligate
toy cross-feeding relationship. SMETANA shares a solver across organisms in
its MUS calculation, while REFRAMED retains minimal-medium enumeration cuts.
Those cuts can constrain a donor when the next organism is evaluated.

`smetana_safe.py` runs the installed SMETANA CLI with a fresh solver for each
organism's MUS minimal-medium calculation. It changes no installed files and
leaves SCS, MIP, and MRO unchanged. This version-specific correction is recorded
in the analysis settings. Always use the supplied pipeline rather than calling
the unwrapped `smetana` command for this environment.

```bash
python smoke_test.py scip
python -m unittest discover -s tests -v
```

The smoke test uses actual SBML models and the actual SCIP/SMETANA code. It
checks monoculture growth, then recovery of the expected formate transfer,
including both donor/receiver input orders. These artificial networks are
software checks and are not biological reconstructions of the LT bins.

The package was tested with those software versions on Linux. It has not been
run on TheBelly or on your actual bin sequences in this session. The production
CarveMe reconstruction requires the input proteins and DIAMOND installation.

## Software sources and methods references

- [CarveMe installation](https://carveme.readthedocs.io/en/latest/installation.html)
  and [usage](https://carveme.readthedocs.io/en/latest/usage.html).
- [SMETANA usage and scores](https://smetana.readthedocs.io/en/latest/usage.html).
- [SMETANA MUS calculation](https://github.com/cdanielmachado/smetana/blob/master/smetana/smetana.py)
  and [REFRAMED minimal-medium enumeration](https://github.com/cdanielmachado/reframed/blob/master/reframed/cobra/medium.py).
- Machado D, Andrejev S, Tramontano M, Patil KR. 2018. Fast automated
  reconstruction of genome-scale metabolic models for microbial species and
  communities. *Nucleic Acids Research* 46:7542–7553.
  [doi:10.1093/nar/gky537](https://doi.org/10.1093/nar/gky537).
- Zelezniak A, Andrejev S, Ponomarova O, Mende DR, Bork P, Patil KR. 2015.
  Metabolic dependencies drive species co-occurrence in diverse microbial
  communities. *PNAS* 112:6449–6454.
  [doi:10.1073/pnas.1421834112](https://doi.org/10.1073/pnas.1421834112).
