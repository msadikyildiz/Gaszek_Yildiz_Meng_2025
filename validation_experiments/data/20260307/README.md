# 20260307 — Experiment Notes

## From Deniz (email)

All 10 plates terminated prematurely but have been read for at least 12 hours.
Main focus: growth up to 12 hours — measurements should be usable as-is.

Discard first 3 OD measurements (plates adjusting to higher temperature and humidity).

Drug: AZT 2916 ug/mL, dilution fold 3.

### Plate layout

| Plate | Variant 1 | Variant 2 |
|-------|-----------|-----------|
| 1     | a2        | a4        |
| 2     | a5        | a6        |
| 3     | a7        | a8        |
| 4     | a9        | a11       |
| 5     | a12       | b1        |
| 6     | b2        | col1      |
| 7     | col2      | col3      |
| 8     | col4      | v1.1      |
| 9     | WT        | Dead (DD) |
| 10    | v2.1      | v3.1      |

### Colony sources
- CML pool: a2, a4, a5, a6, a7, a8, a9, a11, a12
- AZT 4 pool: b1, b2
- Previous (Sanger confirmed): col1-col4
- Known designed variants: v1.1, v2.1, v3.1, WT, DD

## Data quality
- Experiment ran ~13.3h before premature termination
- Valid for 12h analysis; 24h timepoint NOT available
- xlsx has trailing padding rows (time=0, OD=null) — pipeline auto-detects valid range
