# Data provenance and licence

**No data files are committed to this repository.** `make data` downloads them from
the upstream source at a pinned commit. This keeps the question of redistribution
out of this repo entirely.

## What is downloaded

Nine SDF files, one per human cytochrome P450 isoform, from:

- Repository: <https://github.com/dmis-lab/ATTNSOM>
- Path: `cyp_dataset/{1A2,2A6,2B6,2C8,2C9,2C19,2D6,2E1,3A4}.sdf`
- Pinned commit: `0faa38bbd1665262efb3642a2a3ae0f83c928a49` (2026-08-06)
- Repository licence: MIT (`Copyright (c) 2026 DMIS Laboratory - Korea University`)

Contents as downloaded: 2003 molecule–isoform records, 680 unique molecules.
Every record carries `ID`, `PRIMARY_SOM` and `Citation`; 606 carry
`SECONDARY_SOM`, 199 carry `TERTIARY_SOM`, 979 carry a `DOI`. SOM fields hold
1-indexed heavy-atom numbers. Structures are heavy-atom only (hydrogens implicit)
with 3D coordinates from MOE.

Per-isoform record counts, as downloaded:

| 1A2 | 2A6 | 2B6 | 2C8 | 2C9 | 2C19 | 2D6 | 2E1 | 3A4 |
|-----|-----|-----|-----|-----|-----|-----|-----|-----|
| 271 | 105 | 151 | 142 | 226 | 218  | 270 | 145 | 475 |

These match the counts tabulated in the ATTNSOM paper, and the 680-molecule total
matches the set as usually described in the literature.

## Original source of the measurements

The underlying data is the Zaretzki CYP substrate set, originally published as
supporting information to:

- Zaretzki, J.; Bergeron, C.; Rydberg, P.; Huang, T.-W.; Bennett, K. P.; Breneman,
  C. M. *RS-Predictor: A New Tool for Predicting Sites of Cytochrome P450-Mediated
  Metabolism Applied to CYP 3A4.* J. Chem. Inf. Model. **2012**, 52, 1637–1659.
- Zaretzki, J.; Matlock, M.; Swamidass, S. J. *XenoSite: Accurately Predicting
  CYP-Mediated Sites of Metabolism with Neural Networks.* J. Chem. Inf. Model.
  **2013**, 53, 3373–3383. DOI: 10.1021/ci400518g
- Region-resolution extension: Matlock, M.; Hughes, T.; Swamidass, S. J.
  *Extending P450 site-of-metabolism models with region-resolution data.*
  Bioinformatics **2015**, 31, 1966–1973.

Each record's `Citation` / `DOI` field points to the primary experimental paper
the annotation came from.

## Honest statement of the licence position

This is worth being explicit about, because it is the kind of thing that gets
glossed over:

- The **redistribution** in `dmis-lab/ATTNSOM` is under that repository's MIT
  licence, and that is the licence under which this project obtains the files.
- The **original** measurements were published as ACS journal supporting
  information. ACS supporting information is not released under an open licence
  by default, and I have not found an open-licence release from the original
  authors. So the upstream MIT grant covers that repository's distribution of the
  files; it does not retroactively place the underlying compilation under MIT.
- Several independent academic groups redistribute the same set (for example
  `molinfo-vienna/FAME.AL`, MIT, as `data/zaretzki_preprocessed.sdf`), and the set
  is routinely described in the literature as the largest publicly available
  repository of P450 substrate SOM annotations. Academic reuse with citation is
  clearly established practice.

**What this project therefore does:** it ships code only, downloads the data at
run time, and cites the primary sources. It does not re-host the data. Anyone
intending to redistribute the data itself, or to use it commercially, should go
back to the original publishers rather than relying on the MIT grant of a
downstream repository.
