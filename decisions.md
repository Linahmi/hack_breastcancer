# HERoes — decisions log

One block per teammate. Append only; do not edit anyone else's entry.

## Ghita - structure constraints (1N8Z, 4.5 A)
- RED (conserve, direct HER2 contact), 20 positions: 33,50,52,57,58,59,99,103,104,105,163,165,166,167,185,188,226,227,228,229
- GREEN (CDR, mutable): 31,32,34,35,51,53,54,55,56,60,61,62,63,64,65,66,100,101,102,106,107,108,109,159,160,161,162,164,168,169,186,187,189,190,191,224,225,230,231,232
- tool: structure/contacts.py (numpy only) on PDB 1N8Z; chains classified by sequence; numbering asserted residue by residue vs parent.
- reason: preserve the proven paratope (0.82 nM), let PLM scans mutate only non-contact CDR positions (lower risk to affinity).
- correction appliquee: les eaux cristallographiques de la chaine C etaient comptees comme atomes HER2. Corrige -> RED passe de 21 a 20 (la position 201 n avait que ce contact medie par une eau, 5.45 A au plus proche atome proteique) ; GREEN reste 40, inchange, la whitelist de design n a jamais ete affectee. contacts.py --include-waters reproduit l ancien CSV ; commit f01f504 conserve la version pre-correction.
- litterature: 4 hotspots d alanine scan (Kelley 1993) tombent tous sur des positions RED ; la mutation ameliorante publiee D->W (Gerstner 2002, ~3x) tombe sur la position GREEN 102.

## Aida - cofold, gate thresholds, and pool composition

- **Construct.** `inputs/her2_541_620_from_1n8z.pdb` rejected: it is named
  541-620 but contains PDB 541-607 (UniProt 563-629), with a 10-residue
  internal gap at UniProt 603-612 (disordered in 1N8Z chain C) and Cys626
  unpaired because Cys634 falls outside the window. Replaced by
  `src/aida/her2_domainIV.fasta`, UniProt P04626 511-652: 10 intact
  disulfides, zero free cysteines, all 19 Fab-contact residues
  (UniProt 579-627) present with 68/25 aa margin. The C-terminus must
  reach 642 because the 630-642 disulfide closes there, which is why the
  common "domain IV ends at N629" definition was not used.

- **Gate thresholds: KEPT at 0.50 / 0.30. No refit.** Every sample must
  dock at ipSAE >= 0.50 and sample-spread must be <= 0.30. Hedi's
  reasoning, recorded as given: tightening the threshold now that we can
  see where the designs fall would be fitting the threshold to the answer
  we want. The band is empty (closest PASS 0.2608, closest REJECT
  0.3530) and the criterion still separates - it is simply tighter than
  it was on the controls. We state the margin explicitly rather than
  hide it: **the closest passing design sits 0.09 from a rejection.**

- **The threshold question is settled, not worked around.** The only
  borderline point was `lit_hotspot_excluded_01` at spread 0.2608, and
  that design is mislabelled - it carries D102S and so mutates the
  published hotspot its name says it excludes. It is excluded from the
  final pool. No threshold was moved to accommodate anything.

- **The median never gates.** It ranks designs that already passed.
  Demonstrated on a real design: `w01_consensus_green_01` has median
  ipSAE 0.6553 - mid-pack - with a weakest sample of 0.0180 and only 3/5
  docked. A median-only pipeline ships it.

- **`contact_test` keeps all 3 slots despite being the weakest arm.**
  50% pass rate (5/10), and its median design is rejected at 0.3479
  spread. This is a deliberate diversity choice, not an oversight. We
  take its three tightest designs - 0.0577, 0.0674, 0.0873 - because
  dropping a whole strategy from the pool concentrates risk on the arms
  that happened to score well. The arm's weakness is stated, not hidden.

- **Literature arm split by hypothesis.** The blended 60% figure mixed
  two questions and must not be pitched as "published mutations transfer
  poorly". D102W (published mutation): 4/6 PASS, median spread 0.1792.
  Hotspot-excluded (consensus designs with a constraint): 2/4 PASS,
  median spread 0.4390, and effectively n=3 once the mislabelled design
  is removed. Both literature slots test D102W, one coherent hypothesis.
  Against D102S at the same position: 4/6 vs 8/11, Fisher p = 1.000 -
  indistinguishable, n too small to call a winner.

- **esm2_only: absence of a penalty, not superiority.** 12/12 PASS
  across 5-9 mutations within a fixed position set. Fisher vs
  consensus_green 8/10 gives p = 0.195 - not significant. The
  low-perturbation confound is ruled out by the dose-response: 5-mutation
  designs have median spread 0.1020, 7-9 mutation designs 0.0619
  (tighter, not looser), so position choice rather than perturbation size
  is what mattered.

- **Reverse-order control FAILED, reported as a negative result.** Gap
  parent vs reverse 0.0236 against a noise floor of 0.0603-0.0791
  depending on definition - between 0.30x and 0.39x the noise. ipSAE
  cannot resolve the known 2.2x KD difference (0.82 vs 1.83 nM). Checked
  against ipTM (+0.0129 vs 0.0377) and ipSAE asym (+0.0372 vs 0.1081):
  same answer. We do not claim our rankings resolve two-fold affinity.

## Hedi - scoring, controls, and the freeze

### Approche
- Pas de de novo (BindCraft, RFdiffusion, BoltzGen) : on optimise un binder
  qui marche déjà, pas on en invente un.
- Exploration par modèles de séquence (ESM2, AbLang2), cofolding réservé à
  la shortlist. Boltz-2 ne sert qu'à prédire, pas à designer.
- HER2 tronqué au domaine IV (142 aa) : plus rapide, et l'ipSAE (d0res) est
  invariant au trim, contrairement à l'ipTM dont le d0 dépend de la longueur
  des chaînes.

### Contraintes posées
- Linker et cystéines 22/96/158/223 intouchés. Les deux paires tiennent les
  domaines ; une cystéine non appariée = mauvais repliement ou agrégation.
- Numérotation 1-based sur le parent 242 aa : VH 1-120, linker 121-135,
  VL 136-242. Calculée depuis parent.py, jamais codée en dur.

### Le scoring
- ipSAE réimplémenté depuis la référence Dunbrack plutôt que réutilisé tel
  quel. Validé par test de parité : 10170 valeurs synthétiques + un cofold
  réel, 0 écart. Testé par mutation : 9 bugs introduits volontairement,
  9 détectés. Un dixième ne l'a pas été, à juste titre : la variante était
  sémantiquement équivalente.
- Un bug que seule la donnée réelle a révélé : le compte de résidus
  d'interface doit prendre le max des deux directions, pas la direction
  aller. Nos fixtures synthétiques étaient symétriques, donc le test ne
  pouvait pas le voir. Corrigé par une fixture asymétrique.
- Identification des chaînes par la structure (motif du linker), pas par
  pair_chains_iptm dont le mapping vers les IDs YAML est indocumenté.
- 5 diffusion samples par design, agrégation par médiane. La médiane encaisse
  un tirage aberrant ; prendre le meilleur des 5 reviendrait à sélectionner
  le bruit. L'écart entre tirages est conservé comme signal de fiabilité.
  **Amendé après les contrôles — voir la section en fin de document.**

### Les contrôles
- parent_wt (positif), reverse_vl_vh (sonde de topologie), cdrh3_scrambled
  (négatif, composition en acides aminés conservée).
- Lancés à conditions strictement identiques : comparer un scramblé sans MSA
  à un parent avec MSA ne dirait rien.
- Piège évité : dans l'ordre inverse, les cystéines passent à 23/88/144/218.
  Recopier les contraintes du parent aurait créé des ponts sur les mauvais
  résidus, sans que ça se voie.

### Interprétation
- L'ipSAE est un filtre, pas un classement. On écarte les designs cassés,
  puis on choisit 15 designs divers.
- L'ordre inverse ne diffère que d'un facteur 2,2 en KD (0,82 vs 1,83 nM),
  soit ~0,5 kcal/mol — sous la résolution d'un prédicteur de structure.
  Rapporté, non interprété.
- Le vrai contrôle négatif est le CDR-H3 scramblé.

### Gel
16h30. Dernière séquence à 16h30.

### Amendement après les contrôles — agrégation

Les neuf runs de contrôle (3 contrôles × 3 seeds × 5 samples, réglages
identiques) invalident l'agrégation par médiane seule.

`cdrh3_scrambled` n'est pas décalé, il est **bimodal** : un sample soit
s'arrime vers 0,81, soit s'effondre sous 0,07. Aucun des 45 samples ne tombe
entre 0,07 et 0,54. Les effondrements sont réels, pas un artefact de scoring :
l'ipTM_d0chn suit, n0res passe de 242 à 0, l'interface HER2 de 142 à 0, alors
que le pLDDT moyen reste à 82-85 — les deux chaînes se replient, seule
l'interface échoue.

| contrôle | arrimés | médiane | écart intra-seed |
|---|---|---|---|
| parent_wt | 15/15 | 0,9078 | 0,029 – 0,130 |
| reverse_vl_vh | 15/15 | 0,8842 | 0,061 – 0,121 |
| cdrh3_scrambled | 12/15 | 0,8002 | 0,348 – 0,898 |

Chiffres sur la séquence effectivement repliée (97-109
SRWGGDGFYAMDY → RMADDWGFYGSYG, aida 5327aaa). Vérifiés dans les 15 CIF
déposés et les 3 `processed/structures/*.npz`, pas dans `run_meta.json`.
Identiques, au chiffre près, à la dérivation indépendante d'Aida
(`src/aida/controls_ipsae.txt`).

Trois conséquences :

- La médiane d'une variable bimodale rapporte le mode qui a gagné les cinq
  tirages. Les médianes par seed du scramblé vont de 0,8773 à 0,5416, soit
  un empan de 0,3357 pour **une seule séquence à réglages fixes**, contre un
  écart parent–scramblé de 0,1076. Un design jugé sur une seule seed peut
  être retenu par pur tirage.
- La magnitude d'ipSAE ne sépare pas le négatif : quand le scramblé s'arrime,
  il marque 0,8113 contre 0,9050 pour le parent.
- Ce qui sépare, c'est le **taux d'arrimage** (15/15 contre 12/15) et
  l'**écart intra-seed**, disjoints : parent et reverse à 0,1301 au plus,
  scramblé à 0,3481 au moins, facteur 2,7.

**Nouvelle règle.** On filtre sur le taux d'arrimage et l'écart, pas sur la
médiane :

1. tous les samples doivent s'arrimer (ipSAE ≥ 0,5) ;
2. écart intra-run ≤ 0,30 ;
3. la médiane n'intervient qu'ensuite, et seulement entre designs qui
   passent — un petit écart entre deux designs retenus ne devient pas pour
   autant une preuve d'affinité.

Seuils calibrés sur les contrôles : les 6 bons runs passent, les 3 scramblés
échouent. Ils reposent sur n=3 par contrôle et restent provisoires : rien
n'établit le taux d'effondrement réel du parent au-delà de « 0 sur 15 ».

**La marge sur l'écart est mince, et nous ne refittons pas.**

Le seuil est à 0,30. De part et d'autre :

| | écart intra-run | marge au seuil |
|---|---|---|
| contrôle négatif le moins mauvais (`cdrh3_scrambled_seed1`) | 0,3481 | +0,0481 |
| design retenu le plus serré | 0,2608 | −0,0392 |

Le seuil ne tient donc plus qu'à 0,048 du côté négatif et 0,039 du côté
positif. `seed1` arrime 5/5 et n'échoue que sur l'écart : le critère
d'arrimage ne porte plus seul, l'écart est devenu décisif.

**Nous ne recalibrons pas.** Refitter 0,30 maintenant que la distribution est
connue, c'est ajuster le seuil sur la réponse : le seuil cesserait d'être une
prédiction pour devenir une description, et il n'aurait plus aucun pouvoir de
rejet. Le seuil reste à 0,30, tel qu'il a été fixé avant d'avoir vu ces
designs, et la marge est rapportée telle quelle. Si elle est jugée trop
mince, le remède est **plus de contrôles**, pas un seuil déplacé.

Conséquence à assumer : avec 0,039 de marge, un design situé près de la
bordure n'est pas départagé par ce seuil. Les designs concernés doivent être
traités comme indéterminés, pas comme retenus.

L'intuition de garder l'écart comme signal de fiabilité est confirmée. C'est
le meilleur discriminateur des trois ; il devient une barrière, pas un
départage.

#### Deux scramblés, de sévérité différente

Il existe deux contrôles négatifs CDR-H3, et ils ne se valent pas :

| scramblé | fenêtre | permutation | arrimés | médiane | écart min |
|---|---|---|---|---|---|
| A (Aida) | 99-109 | WGGDGFYAMDY → GFDMAWYGYDG | 7/15 | 0,2754 | 0,832 |
| B (Hedi) | 97-109 | SRWGGDGFYAMDY → RMADDWGFYGSYG | 12/15 | 0,8002 | 0,348 |

B est le moins destructeur des deux, malgré une fenêtre plus large : il
conserve `GFYG` en place et l'interface se reforme la plupart du temps.

**Le seuil de 0,30 a été calibré sur A**, le plus sévère, dont l'écart
minimal était à 0,832 — d'où la marge confortable de +0,53 à l'époque. Il est
aujourd'hui évalué contre B, à 0,348, d'où les +0,048 ci-dessus. Le seuil n'a
pas bougé ; c'est le contrôle qui a changé, parce que les répertoires
`cdrh3_scrambled_*` ont été re-roulés depuis le YAML de B (aida 5327aaa).

Seul B est actuellement sur disque : les fichiers de A ont été remplacés par
le re-run, et ses chiffres ci-dessus proviennent du scoring effectué sur le
snapshot dac6c33. Relancer A comme second négatif rétablirait une marge
large, et deux négatifs de sévérité inégale renseignent plus qu'un seul.

Note méthodologique sur le plancher de bruit : Aida rapporte l'écart
parent–scramblé contre l'écart inter-seed **du parent** (0,0603), ce qui
donne 1,78× et « résolu ». Nous rapportons aussi le même écart contre celui
du scramblé (0,3357), ce qui donne 0,3×. Les deux sont justes et répondent à
deux questions différentes : « l'écart dépasse-t-il le bruit d'un design
sain ? » (oui) et « une seule seed suffit-elle à trancher ? » (non). Les
deux chiffres bruts sont identiques entre nos deux dérivations.

## Aida — POST-HOC CHECK (run after selection; did not change it)

Labelled post-hoc deliberately. The 15 were already chosen and committed
when these numbers were produced. `decisions.md` commits us, before any
design score existed, to selection by diversity across arms plus sample
spread within them. Nothing below was used to re-select, and no design
was swapped. The question asked was only whether a different signal
agrees with the choice.

- **pDockQ2 and LIS recovered without new compute.** Hedi's `ipsae.py`
  computes them (they are in its `COLUMNS`), but `score_designs.py` read
  only `ipsae` and discarded the rest. Re-reading the 74 existing run
  directories recovered them. No cofolding.

- **Correlation with ipSAE median, across the 43 passing designs:**
  pDockQ2 Pearson +0.649, Spearman +0.506. LIS Pearson +0.883,
  Spearman +0.921. pDockQ Pearson +0.493, Spearman +0.273.

- **A prediction recorded before these numbers existed was WRONG.** Aida
  predicted pDockQ2 would correlate at 0.8 or above with ipSAE, on the
  reasoning that both derive from the same PAE matrix. It came back at
  0.649 Pearson / 0.506 Spearman. pDockQ2 is more independent of ipSAE
  than expected. LIS, at 0.883/0.921, is the metric that behaved as
  predicted — and is therefore the one carrying the least new
  information.

- **The chosen 15 are spread across the range, not clustered at the top.**
  On pDockQ2 their ranks among the 43 are 2, 3, 5, 6, 10, 12, 14, 15, 17,
  21, 29, 30, 33, 34, 36 — median rank 15 of 43, with 8 in the top 15 and
  10 in the top half. On LIS: 1, 2, 3, 5, 7, 8, 12, 15, 17, 21, 25, 26,
  29, 32, 38, also median rank 15. Chosen mean pDockQ2 0.5187 against
  0.4840 for the 28 not picked; the single best pDockQ2 in the pool
  (0.6758) is NOT in our 15.

  So an orthogonal-ish signal **mildly agrees and does not endorse**. Our
  pool is better than average on pDockQ2 but we did not pick its top. That
  is the expected consequence of selecting on consistency and arm
  diversity rather than on score magnitude, and it is the honest answer:
  had we optimised pDockQ2 we would have chosen differently.

- **Mean pLDDT, parsed from data already in `scored.csv`.** Parent mean
  87.21 over 15 samples (range 85.7–88.8, per-seed means 87.40/87.42/86.80,
  so seed-level spread 0.62). 13 of our 15 sit below the parent mean, span
  85.72–87.34, largest deficit −1.49 (`esm2_only_07`). All 15 fall inside
  the parent's own per-sample range, but the largest deficits exceed the
  parent's seed-to-seed spread of 0.62, so they are not purely noise.
  Flagged as a mild expression risk worth knowing, since expression is
  published on ProteinBase. Not acted on: a 1.5-point pLDDT difference is
  not a basis for swapping a design, and we had committed to the
  criterion in advance.
