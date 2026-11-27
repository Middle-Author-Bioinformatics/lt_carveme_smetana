#!/usr/bin/env python3
"""Reconstruct major LT bins with CarveMe; run SMETANA within each culture.

Input is per-bin PROTEIN FASTA (.faa/.faa.gz), not contig DNA or annotation TSV.
See README.md for installation, medium definitions, and interpretation.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import importlib.metadata
import itertools
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

VERSION = '1.0.0'
HERE = Path(__file__).resolve().parent
EXPECTED = {'LT1': {'Geobacter'}, 'LT2': {'Geobacter'},
            'LT3': {'Geobacter', 'Azospira'},
            'LT4': {'Geobacter', 'Azospira', 'Cupidesulfovibrio'},
            'LT5': {'Geobacter', 'Azospira', 'Cupidesulfovibrio'}}


def write_tsv(path, fields, rows):
    with Path(path).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with Path(path).open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def inventory(args):
    rows = []
    if args.manifest:
        manifest = args.manifest.resolve()
        rows = read_tsv(manifest)
        required = {'sample', 'taxon', 'fasta'}
        if not rows or not required.issubset(rows[0]):
            raise ValueError('Manifest needs a header: sample\ttaxon\tfasta')
        for row in rows:
            row['fasta'] = str((manifest.parent / row['fasta']).resolve())
    else:
        if not args.proteins.is_dir():
            raise ValueError(f'Protein directory does not exist: {args.proteins}')
        files = sorted(p for p in args.proteins.iterdir() if p.is_file() and p.name.lower().endswith(('.faa', '.faa.gz')))
        for path in files:
            match = re.search(r'(?<![A-Za-z0-9])(LT[1-5])(?![0-9])', path.name, re.I)
            taxon = next((t for t in ['Nitratidesulfovibrio', 'Cupidesulfovibrio', 'Geobacter', 'Azospira']
                          if t.lower() in path.name.lower()), None)
            if not match or not taxon:
                print(f'Excluded from major-bin analysis: {path.name}', flush=True)
                continue
            rows.append({'sample': match[1].upper(), 'taxon': taxon, 'fasta': str(path.resolve())})
    seen = set()
    for row in rows:
        sample = row['sample'].strip().upper()
        taxon = row['taxon'].strip()
        canonical = {'nitratidesulfovibrio': 'Cupidesulfovibrio',
                     'cupidesulfovibrio': 'Cupidesulfovibrio', 'geobacter': 'Geobacter', 'azospira': 'Azospira'}
        taxon = canonical.get(taxon.lower(), taxon)
        if sample not in EXPECTED or taxon not in EXPECTED[sample]:
            raise ValueError(f'Unexpected major bin {sample}/{taxon}; check the manifest.')
        if (sample, taxon) in seen:
            raise ValueError(f'Multiple files for {sample}/{taxon}; use a manifest to choose one.')
        if not Path(row['fasta']).is_file():
            raise ValueError(f'Missing FASTA: {row["fasta"]}')
        seen.add((sample, taxon))
        row.update(sample=sample, taxon=taxon, model_id=f'{sample}_{taxon}')
    missing = sorted((s, t) for s, taxa in EXPECTED.items() for t in taxa if (s, t) not in seen)
    if missing:
        raise ValueError('Missing expected major bins: ' + ', '.join('/'.join(x) for x in missing))
    return sorted(rows, key=lambda r: r['model_id'])


def normalize_proteins(source, target, mapping):
    """Use stable safe IDs for CarveMe and preserve the original FASTA headers."""
    opener = gzip.open if str(source).lower().endswith('.gz') else open
    records, header, sequence = [], None, []
    def add_record():
        if header is None:
            return
        seq = ''.join(sequence).upper().rstrip('*')
        if not seq or set(seq) - set('ACDEFGHIKLMNPQRSTVWYBXZJUO'):
            raise ValueError(f'Invalid protein sequence in {source}: {header}')
        if len(seq) > 60 and set(seq) <= set('ACGTN'):
            raise ValueError(f'Looks like DNA, not protein: {source}: {header}')
        records.append((header, seq))
    with opener(source, 'rt') as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                add_record()
                header, sequence = line[1:], []
            elif header is None:
                raise ValueError(f'Expected a FASTA header in {source}')
            else:
                sequence.append(line)
        add_record()
    if not records:
        raise ValueError(f'No proteins in {source}')
    with target.open('w') as handle:
        for number, (_, seq) in enumerate(records, 1):
            handle.write(f'>gene_{number:07d}\n{seq}\n')
    write_tsv(mapping, ['model_gene_id', 'original_header'],
              [{'model_gene_id': f'gene_{i:07d}', 'original_header': h} for i, (h, _) in enumerate(records, 1)])
    return len(records)


def run_command(command, log):
    with Path(log).open('w') as handle:
        handle.write(shlex.join(map(str, command)) + '\n\n')
        handle.flush()
        result = subprocess.run(list(map(str, command)), stdout=handle, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f'Command failed ({result.returncode}); see {log}')


def load_model(path):
    from reframed import load_cbmodel
    model = load_cbmodel(str(path), flavor='fbc2')
    if not model.get_objective() or not model.biomass_reaction:
        raise ValueError(f'Model lacks a biomass objective: {path}')
    # Match SMETANA interface.build_cache, which disables ATP maintenance.
    if 'R_ATPM' in model.reactions:
        model.reactions['R_ATPM'].lb = 0
    return model


def reconstruct(row, args, versions):
    name, root = row['model_id'], args.outdir
    xml = root/'models'/f'{name}.xml'
    metadata = xml.with_suffix('.json')
    signature = {'pipeline': VERSION, 'input_sha256': sha(row['fasta']), 'solver': args.solver,
                 'versions': versions, 'universe': 'bacteria', 'extra_medium_gapfill': False}
    if xml.exists() or metadata.exists():
        old = json.loads(metadata.read_text()) if metadata.exists() else {}
        if args.resume and old.get('signature') == signature and xml.exists() and old.get('model_sha256') == sha(xml):
            load_model(xml)
            print(f'Reusing checked model: {name}', flush=True)
            return xml
        raise ValueError(f'{xml} already exists or its settings changed. Use --resume for matching models or a new output directory.')
    staged = root/'inputs'/f'{name}.faa'
    count = normalize_proteins(row['fasta'], staged, root/'inputs'/f'{name}.gene_ids.tsv')
    # A separate directory retains CarveMe/DIAMOND intermediates and partial outputs.
    work = root/'work'/name
    work.mkdir(parents=True, exist_ok=True)
    temporary = work/f'{name}.xml'
    temporary.unlink(missing_ok=True)
    command = ['carve', staged, '-o', temporary, '--fbc2', '--solver', args.solver,
               '--diamond-args', f'--more-sensitive --top 10 --quiet --threads {args.threads}', '-v']
    print(f'Carving {name} ({count} proteins)', flush=True)
    run_command(command, root/'logs'/f'{name}.carveme.log')
    if not temporary.is_file() or temporary.stat().st_size == 0:
        raise RuntimeError(f'CarveMe produced no SBML for {name}; inspect its log.')
    load_model(temporary)
    temporary.replace(xml)
    metadata.write_text(json.dumps({'signature': signature, 'model_sha256': sha(xml), 'n_proteins': count,
                                    'command': list(map(str, command))}, indent=2)+'\n')
    return xml


def communities(rows, pairs):
    groups = {}
    for sample in EXPECTED:
        ids = [r['model_id'] for r in rows if r['sample'] == sample]
        if len(ids) < 2:
            continue
        groups[sample] = ids
        if pairs and len(ids) > 2:
            for pair in itertools.combinations(ids, 2):
                groups[sample+'__'+ '__'.join(p.removeprefix(sample+'_') for p in pair)] = list(pair)
    return groups


def media_settings(args):
    if bool(args.media_db) != bool(args.media):
        raise ValueError('--media-db and --media must be supplied together.')
    if args.exploratory and args.media:
        raise ValueError('Choose --exploratory or a supplied medium.')
    if not args.models_only and not args.plan_only and not (args.exploratory or args.media):
        raise ValueError('Supply --media-db / --media, or explicitly select --exploratory. See README.md.')
    result = {}
    if args.media:
        for row in read_tsv(args.media_db):
            if not {'medium', 'compound'}.issubset(row):
                raise ValueError('Media TSV needs a header: medium\tcompound')
            name, compound = row['medium'].strip(), row['compound'].strip()
            if not name or not re.fullmatch(r'[A-Za-z0-9_]+', compound):
                raise ValueError('Use bare BiGG compound IDs in medium: ac, fum, h2o, etc.')
            result.setdefault(name, set()).add(compound)
        selected = args.media.split(',')
        if any(not x or x not in result for x in selected):
            raise ValueError('Every requested medium must occur in --media-db.')
        if any('o2' in result[x] for x in selected):
            raise ValueError('A selected anaerobic medium contains o2; remove its row.')
        result = {x: sorted(result[x]) for x in selected}
    return result


def analyze(rows, paths, args, media):
    from reframed import Environment, FBA
    from reframed.solvers.solution import Status
    from smetana.interface import define_environment
    from smetana.legacy import Community
    analysis = args.outdir/'analysis'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    analysis.mkdir(parents=True)
    print(f'Analysis: {analysis}', flush=True)
    (analysis/'settings.json').write_text(json.dumps({'argv': sys.argv, 'pipeline': VERSION,
        'exploratory': args.exploratory, 'media': media, 'media_sha256': sha(args.media_db) if args.media_db else None,
        'MUS_solver_isolation': True, 'anaerobic': True, 'SMETANA_min_growth': 0.1,
        'model_sha256': {p.stem: sha(p) for p in paths}}, indent=2)+'\n')
    if media:
        write_tsv(analysis/'media.tsv', ['medium', 'compound'],
                  [{'medium': m, 'compound': c} for m, compounds in media.items() for c in compounds])
    models = {p.stem: load_model(p) for p in paths}
    exchange_rows, growth_rows = [], []
    for name, model in models.items():
        for reaction_id in model.get_exchange_reactions():
            reaction = model.reactions[reaction_id]
            for met in reaction.stoichiometry:
                exchange_rows.append({'model_id': name, 'reaction': reaction_id, 'metabolite': met,
                                      'name': model.metabolites[met].name})
        conditions = media or {'complete_anaerobic_exploratory': None}
        for medium, compounds in conditions.items():
            copy = model.copy()
            env = Environment.from_compounds(compounds) if compounds is not None else Environment.complete(copy)
            if compounds is None and 'R_EX_o2_e' in copy.reactions:
                env['R_EX_o2_e'] = (0, math.inf)
            env.apply(copy, inplace=True, warning=False)
            solution = FBA(copy, get_values=False)
            growth_rows.append({'model_id': name, 'medium': medium, 'status': str(solution.status),
                                'maximum_biomass': solution.fobj,
                                'reaches_smetana_threshold': solution.status == Status.OPTIMAL and solution.fobj >= 0.1})
    write_tsv(analysis/'model_exchanges.tsv', ['model_id', 'reaction', 'metabolite', 'name'], exchange_rows)
    write_tsv(analysis/'monoculture_growth.tsv', ['model_id', 'medium', 'status', 'maximum_biomass', 'reaches_smetana_threshold'], growth_rows)
    # Missing requested substrates are visible even when FBA cannot grow.
    absent = []
    for name, model in models.items():
        available = {m.removeprefix('M_').removesuffix('_e') for r in model.get_exchange_reactions()
                     for m in model.reactions[r].stoichiometry if m.endswith('_e')}
        for medium, compounds in media.items():
            absent.extend({'model_id': name, 'medium': medium, 'compound': c} for c in compounds if c not in available)
    write_tsv(analysis/'medium_compounds_without_exchange.tsv', ['model_id', 'medium', 'compound'], absent)
    combined = {'global': [], 'detailed': []}
    status_rows = []
    for group, members in communities(rows, args.pairs).items():
        group_dir = analysis/group
        group_dir.mkdir()
        membership = group_dir/'community.tsv'
        membership.write_text(''.join(f'{group}\t{x}\n' for x in members))
        for mode in ('global', 'detailed'):
            for medium in media or {None: None}:
                # Preflight the SAME environment used by this SMETANA mode.
                community = Community(group, [models[x].copy() for x in members], copy_models=False)
                medium_name, env = define_environment(medium, media, community, mode, False, False, False, False)
                env.apply(community.merged, inplace=True, warning=False)
                biomass = list(community.organisms_biomass_reactions.values())
                solution = FBA(community.merged, objective={x: 1 for x in biomass},
                               constraints={x: (0.1, math.inf) for x in biomass}, get_values=False)
                condition = (medium if medium else 'exploratory')
                # Paths use an index to keep arbitrary medium labels out of filenames.
                prefix = group_dir/f'{mode}_{len(status_rows):03d}'
                record = {'community': group, 'mode': mode, 'medium': medium_name,
                          'joint_growth_status': str(solution.status), 'result': 'skipped_no_joint_growth',
                          'output': '', 'log': ''}
                if solution.status == Status.OPTIMAL:
                    command = [sys.executable, HERE/'smetana_safe.py', *[args.outdir/'models'/f'{x}.xml' for x in members],
                               '-c', membership, '--flavor', 'fbc2', '--solver', args.solver, '--anaerobic',
                               '-g' if mode == 'global' else '-d', '-o', prefix, '-v', '-z']
                    if medium:
                        command += ['-m', medium, '--mediadb', analysis/'media.tsv']
                    print(f'SMETANA {group}: {mode}, {condition}', flush=True)
                    log = Path(str(prefix)+'.log')
                    run_command(command, log)
                    output = Path(str(prefix)+f'_{mode}.tsv')
                    result = read_tsv(output)
                    combined[mode].extend(result)
                    # Retain warnings even when SMETANA exits successfully.
                    text = log.read_text()
                    has_warning = any(x in text.lower() for x in ('warning', 'failed to', 'no solution'))
                    missing_values = any(r.get(k, '') in ('', 'nan', 'None') for r in result
                                         for k in (['mip', 'mro'] if mode == 'global' else ['scs', 'mus', 'mps', 'smetana']))
                    record.update(result='review_log' if has_warning or missing_values else ('ok' if result else 'empty_review_log'),
                                  output=str(output.relative_to(analysis)), log=str(log.relative_to(analysis)))
                else:
                    print(f'Skipping {group}/{mode}/{medium_name}: joint growth infeasible; see run_status.tsv', flush=True)
                status_rows.append(record)
                write_tsv(analysis/'run_status.tsv', list(record), status_rows)
    fields = {'global': ['community', 'medium', 'size', 'mip', 'mro'],
              'detailed': ['community', 'medium', 'receiver', 'donor', 'compound', 'scs', 'mus', 'mps', 'smetana']}
    for mode in combined:
        write_tsv(analysis/f'all_{mode}.tsv', fields[mode], combined[mode])
    positive = sorted((r for r in combined['detailed'] if float(r['smetana']) > 0), key=lambda r: -float(r['smetana']))
    write_tsv(analysis/'positive_exchanges.tsv', fields['detailed'], positive)
    print(f'Done. Candidate exchanges: {analysis}/positive_exchanges.tsv\nCheck model growth and run_status.tsv before interpreting scores.', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--proteins', type=Path, help='Folder containing the ten major-bin LT*.faa[.gz] files')
    source.add_argument('--manifest', type=Path, help='TSV with sample, taxon, fasta columns; paths relative to manifest')
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--jobs', type=int, default=4, help='Concurrent CarveMe jobs (default 4)')
    parser.add_argument('--threads', type=int, default=8, help='DIAMOND threads per job (default 8)')
    parser.add_argument('--solver', choices=['scip', 'gurobi', 'cplex'], default='scip')
    parser.add_argument('--media-db', type=Path, help='TSV with medium, compound columns; bare BiGG IDs')
    parser.add_argument('--media', help='Comma-separated medium names; all models use each selected medium')
    parser.add_argument('--exploratory', action='store_true', help='Use SMETANA model-derived anaerobic defaults; not an assay simulation')
    parser.add_argument('--pairs', action='store_true', help='Also run every pair within LT4 and LT5')
    parser.add_argument('--resume', action='store_true', help='Reuse matching reconstructed models; analysis gets a new dated directory')
    parser.add_argument('--models-only', action='store_true')
    parser.add_argument('--plan-only', action='store_true', help='Validate inputs and save resolved manifest without external programs')
    args = parser.parse_args()
    os.environ.setdefault('OMP_NUM_THREADS', '1')
    os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
    os.environ.setdefault('MKL_NUM_THREADS', '1')
    if args.jobs < 1 or args.threads < 1:
        parser.error('--jobs and --threads must be positive')
    args.outdir = args.outdir.resolve()
    rows = inventory(args)
    media = media_settings(args)
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.outdir/'resolved_manifest.tsv', ['sample', 'taxon', 'model_id', 'fasta'], rows)
    if args.plan_only:
        print(f'Validated {len(rows)} bins. See {args.outdir}/resolved_manifest.tsv')
        return
    for program in ('carve', 'smetana', 'diamond'):
        if not shutil.which(program):
            raise RuntimeError(f'{program} is missing. Activate the installed environment.')
    versions = {p: importlib.metadata.version(p) for p in ['carveme', 'smetana', 'reframed', 'pyscipopt']}
    versions['diamond'] = subprocess.check_output(['diamond', 'version'], text=True).strip()
    if any(versions[p] != v for p, v in [('carveme', '1.6.6'), ('smetana', '1.2.1'), ('reframed', '1.6.0')]):
        raise RuntimeError(f'Use the pinned software environment: {versions}')
    from reframed import set_default_solver
    set_default_solver(args.solver)
    for directory in ('models', 'inputs', 'work', 'logs'):
        (args.outdir/directory).mkdir(exist_ok=True)
    run_command([sys.executable, '-m', 'pip', 'freeze'], args.outdir/'logs'/'pip_freeze.txt')
    run_command(['diamond', 'version'], args.outdir/'logs'/'diamond_version.txt')
    print(f'{args.jobs} model jobs; {args.threads} DIAMOND threads each ({args.jobs*args.threads} total).', flush=True)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        paths = list(pool.map(lambda row: reconstruct(row, args, versions), rows))
    if not args.models_only:
        analyze(rows, paths, args, media)
    else:
        print(f'Done. SBML models: {args.outdir}/models')


if __name__ == '__main__':
    try:
        main()
    except (ValueError, RuntimeError, OSError, ImportError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
