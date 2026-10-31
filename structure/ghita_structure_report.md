# HERoes - Structure & constraints (Ghita)
Source: PDB 1N8Z (HER2 ectodomain + trastuzumab Fab, 2.52 A). Chains auto-detected: A=VL, B=VH, C=HER2. Contact cutoff 4.5 A (any atom).
Numbering verified: every mapped PDB residue matches the parent 1-242 sequence.

NOTE for the team: brief boundaries are off by one. TRUE: VH 1-120, linker 121-135, VL 136-242.

CORRECTION (see structure/README.md): the first version of this report counted crystallographic
waters of chain C as HER2 atoms, which gave 21 RED and a 17-residue core at 4.0 A. Waters are now
excluded. RED is 20 - position 201 (framework) was RED on a single water and has no HER2 protein
atom within 4.5 A (nearest 5.45 A) - and the 4.0 A core is 14 residues. GREEN is unchanged at 40.
Distances for 58 and 163 moved to 4.11 and 4.50 A. Regenerate with:
    python3 structure/contacts.py 1N8Z.cif > structure/red_green.csv

## RED - contacts HER2, DO NOT mutate (20)
-  33 Y  CDR-H1  (3.42 A)
-  50 R  CDR-H2  (2.60 A)
-  52 Y  CDR-H2  (4.39 A)
-  57 Y  CDR-H2  (3.51 A)
-  58 T  CDR-H2  (4.11 A)
-  59 R  CDR-H2  (3.43 A)
-  99 W  CDR-H3  (4.05 A)
- 103 G  CDR-H3  (2.42 A)
- 104 F  CDR-H3  (4.27 A)
- 105 Y  CDR-H3  (3.32 A)
- 163 D  CDR-L1  (4.50 A)
- 165 N  CDR-L1  (2.99 A)
- 166 T  CDR-L1  (3.05 A)
- 167 A  CDR-L1  (3.69 A)
- 185 S  CDR-L2  (3.61 A)
- 188 F  CDR-L2  (3.56 A)
- 226 H  CDR-L3  (2.89 A)
- 227 Y  CDR-L3  (2.96 A)
- 228 T  CDR-L3  (4.05 A)
- 229 T  CDR-L3  (2.82 A)

## GREEN - CDR positions, safe to mutate (40)
31, 32, 34, 35, 51, 53, 54, 55, 56, 60, 61, 62, 63, 64, 65, 66, 100, 101, 102, 106, 107, 108, 109, 159, 160, 161, 162, 164, 168, 169, 186, 187, 189, 190, 191, 224, 225, 230, 231, 232

## By CDR (green = mutable)
- CDR-H1 (31-35): RED [33] | GREEN [31, 32, 34, 35]
- CDR-H2 (50-66): RED [50, 52, 57, 58, 59] | GREEN [51, 53, 54, 55, 56, 60, 61, 62, 63, 64, 65, 66]
- CDR-H3 (99-109): RED [99, 103, 104, 105] | GREEN [100, 101, 102, 106, 107, 108, 109]
- CDR-L1 (159-169): RED [163, 165, 166, 167] | GREEN [159, 160, 161, 162, 164, 168, 169]
- CDR-L2 (185-191): RED [185, 188] | GREEN [186, 187, 189, 190, 191]
- CDR-L3 (224-232): RED [226, 227, 228, 229] | GREEN [224, 225, 230, 231, 232]
