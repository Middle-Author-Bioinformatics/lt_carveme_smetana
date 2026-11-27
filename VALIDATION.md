# Validation

Checked on Linux with Python 3.12, CarveMe 1.6.6, SMETANA 1.2.1,
REFRAMED 1.6.0, and PySCIPOpt 6.3.0.

- Shell syntax check: passed.
- Python compilation: passed.
- Seven unit/integration tests: passed. These cover input selection and the
  revised taxon name, missing/duplicate bin rejection, within-culture grouping,
  protein/DNA validation and gene-ID mapping, explicit anaerobic medium
  selection, detection of an apparently successful CarveMe command that produces
  no SBML, and a real SBML-to-SMETANA analysis with growing and infeasible media.
- The real solver smoke test passes with both donor-first and receiver-first
  community order. The donor grows on the supplied substrate; the receiver
  fails to grow alone; the corrected detailed calculation recovers the expected
  donor-to-receiver formate exchange. Global calculations also execute.
- The same donor-first toy case produces no positive detailed rows using the
  unmodified SMETANA CLI. The supplied wrapper addresses the shared-solver
  enumeration-cut carryover described in README.md.

No actual LT bin files were available for a production reconstruction here.
The conda installer and DIAMOND/CarveMe reconstruction have not been executed
on TheBelly. The installer runs the solver smoke test on the target machine.
