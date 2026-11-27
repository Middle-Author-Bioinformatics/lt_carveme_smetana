#!/usr/bin/env python3
"""SMETANA CLI with isolated MUS enumeration for the pinned software versions.

SMETANA 1.2.1 mu_score shares a solver between species. REFRAMED 1.6.0
minimal_medium leaves enumeration cuts in that solver. Those cuts can prevent
a later species from using an earlier species as a donor. Keep the published
algorithm, but start each MUS minimal-medium calculation with a fresh solver.
This wrapper changes no installed files and does not affect SCS, MIP, or MRO.
"""
import importlib.metadata
import runpy
import shutil


def apply_fix():
    import smetana.interface as interface
    import smetana.smetana as algorithm
    from reframed import solver_instance
    for package, expected in [('smetana', '1.2.1'), ('reframed', '1.6.0')]:
        actual = importlib.metadata.version(package)
        if actual != expected:
            raise RuntimeError(f'This wrapper was checked with {package} {expected}; found {actual}. Use the supplied environment.')
    original_mu = algorithm.mu_score

    def isolated_mu(*args, **kwargs):
        original_minimal = algorithm.minimal_medium

        def fresh_minimal(model, *medium_args, **medium_kwargs):
            medium_kwargs['solver'] = solver_instance(model)
            return original_minimal(model, *medium_args, **medium_kwargs)

        algorithm.minimal_medium = fresh_minimal
        try:
            return original_mu(*args, **kwargs)
        finally:
            algorithm.minimal_medium = original_minimal

    interface.mu_score = isolated_mu


if __name__ == '__main__':
    apply_fix()
    executable = shutil.which('smetana')
    if not executable:
        raise SystemExit('smetana is missing; activate the supplied environment.')
    runpy.run_path(executable, run_name='__main__')
