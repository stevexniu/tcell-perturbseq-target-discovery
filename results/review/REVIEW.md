# Independent review — verdict: RELEASE-READY

A fresh-context reviewer traced every per-gene claim in the reports against
`results/tables/coverage_table.csv` and inspected all three figures.

**7 clean sections** (per-gene claims row-for-row; direction-of-effect logic incl.
CBLB risk-on-LoF; activation-axis weakness disclosed everywhere; 11/15 & 13/15
aggregates reconcile over all 15 candidates; Figs 1-3 correct).

**2 prior fixes confirmed resolved:** Fig-2 label stacks ordered by descending
score (no crossing leaders); CLCC1 (15th candidate) present in both reports and
the coverage table.

**3 open items addressed in the final report:**
- (medium) CLCC1 was described as "among the more T-cell-specific" — corrected:
  its K562 concordance (0.48) is the highest/least-specific of the 8 measured
  candidates.
- (minor) Fig-3 "direction consistent across conditions" — qualified: holds for
  the polarization nominations; CD6 flips sign by 48 h and INTS6/TATDN2 decay to ~0.
- (info) INTS6/SNRPB2 carry registered-but-empty ChEMBL target IDs while called
  genuinely-new — no change needed; the "as a target / no drugs" qualifier is present.

Full machine-readable findings: `review_findings.json`.
