#!/usr/bin/env python3
"""Exercise real SBML, FBA, and SMETANA with an obligate toy cross-feeding pair."""
import csv
from pathlib import Path
import subprocess
import sys
import tempfile


def make_toy(path, donor):
    from reframed import CBModel, CBReaction, Compartment, Metabolite, save_cbmodel
    model = CBModel(path.stem)
    model.add_compartment(Compartment('c'))
    model.add_compartment(Compartment('e', external=True))
    for compound in ('ac', 'for', 'o2', 'h2o'):
        for compartment in ('c', 'e'):
            model.add_metabolite(Metabolite(f'M_{compound}_{compartment}', compartment=compartment))
        model.add_reaction(CBReaction(f'R_EX_{compound}_e', stoichiometry={f'M_{compound}_e': -1}, lb=-10, ub=1000))
        model.add_reaction(CBReaction(f'R_T_{compound}', stoichiometry={f'M_{compound}_e': -1, f'M_{compound}_c': 1}, lb=-1000, ub=1000))
    model.add_metabolite(Metabolite('M_precursor_c', compartment='c'))
    stoich = {'M_ac_c': -1, 'M_precursor_c': 1, 'M_for_c': 1} if donor else {'M_for_c': -1, 'M_precursor_c': 1}
    model.add_reaction(CBReaction('R_SYNTHESIS', stoichiometry=stoich, lb=0, ub=1000))
    model.add_reaction(CBReaction('R_BIOMASS', stoichiometry={'M_precursor_c': -1}, lb=0, ub=1000, objective=1))
    save_cbmodel(model, str(path), flavor='fbc2')


def main():
    solver = sys.argv[1] if len(sys.argv) > 1 else 'scip'
    from reframed import FBA, Environment, load_cbmodel, set_default_solver
    set_default_solver(solver)
    with tempfile.TemporaryDirectory(prefix='lt-smoke-') as tmp:
        root = Path(tmp)
        donor, receiver = root/'donor.xml', root/'receiver.xml'
        make_toy(donor, True)
        make_toy(receiver, False)
        (root/'media.tsv').write_text('medium\tcompound\ntoy\tac\ntoy\th2o\n')
        (root/'communities.tsv').write_text('toy\tdonor\ntoy\treceiver\n')
        for path, expected_growth in [(donor, True), (receiver, False)]:
            model = load_cbmodel(str(path), flavor='fbc2')
            Environment.from_compounds(['ac', 'h2o']).apply(model, inplace=True, warning=False)
            solution = FBA(model)
            assert (solution.fobj is not None and solution.fobj > 1e-6) == expected_growth
        for members in [('donor', 'receiver'), ('receiver', 'donor')]:
            (root/'communities.tsv').write_text(''.join(f'toy\t{x}\n' for x in members))
            for mode in ('-g', '-d'):
                subprocess.run([sys.executable, str(Path(__file__).with_name('smetana_safe.py')), str(donor), str(receiver), '-c', str(root/'communities.tsv'),
                                '-m', 'toy', '--mediadb', str(root/'media.tsv'), '--solver', solver,
                                '--flavor', 'fbc2', '--anaerobic', mode, '-o', str(root/'toy')], check=True)
            with (root/'toy_detailed.tsv').open() as handle:
                rows = list(csv.DictReader(handle, delimiter='\t'))
            assert any(r['donor'] == 'donor' and r['receiver'] == 'receiver'
                       and r['compound'] == 'M_for_e' and float(r['smetana']) > 0 for r in rows), rows
    print(f'PASS: {solver}; donor grows alone, receiver needs donor, SMETANA recovers formate exchange.')


if __name__ == '__main__':
    main()
