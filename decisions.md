# HERoes — decisions log

One block per teammate. Append only; do not edit anyone else's entry.

## Ghita - structure constraints (1N8Z, 4.5 A)
- RED (conserve, direct HER2 contact): 33,50,52,57,58,59,99,103,104,105,163,165,166,167,185,188,201,226,227,228,229
- GREEN (CDR, mutable): 31,32,34,35,51,53,54,55,56,60,61,62,63,64,65,66,100,101,102,106,107,108,109,159,160,161,162,164,168,169,186,187,189,190,191,224,225,230,231,232
- tool: Bio.PDB NeighborSearch on PDB 1N8Z; chains classified by sequence; numbering checked vs parent.
- reason: preserve the proven paratope (0.82 nM), let PLM scans mutate only non-contact CDR positions (lower risk to affinity).
- correction: les eaux cristallographiques de la chaine C etaient comptees comme atomes HER2. Sans elles, RED = 20 (la position 201 perd son seul contact, 5.45 A au plus proche atome proteique) ; GREEN reste 40, inchange. contacts.py --include-waters reproduit le CSV original.
- litterature: 4 hotspots d alanine scan (Kelley 1993) tombent tous sur des positions RED ; la mutation ameliorante publiee D->W (Gerstner 2002, ~3x) tombe sur la position GREEN 102.
