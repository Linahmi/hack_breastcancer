# decisions.md — HERoes

## Approche
- Pas de de novo (BindCraft, RFdiffusion, BoltzGen) : on optimise un binder
  qui marche déjà, pas on en invente un.
- Exploration par modèles de séquence (ESM2, AbLang2), cofolding réservé à
  la shortlist. Boltz-2 ne sert qu'à prédire, pas à designer.
- HER2 tronqué au domaine IV (142 aa) : plus rapide, et l'ipSAE (d0res) est
  invariant au trim, contrairement à l'ipTM dont le d0 dépend de la longueur
  des chaînes.

## Contraintes posées
- Linker et cystéines 22/96/158/223 intouchés. Les deux paires tiennent les
  domaines ; une cystéine non appariée = mauvais repliement ou agrégation.
- Numérotation 1-based sur le parent 242 aa : VH 1-120, linker 121-135,
  VL 136-242. Calculée depuis parent.py, jamais codée en dur.

## Le scoring
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

## Les contrôles
- parent_wt (positif), reverse_vl_vh (sonde de topologie), cdrh3_scrambled
  (négatif, composition en acides aminés conservée).
- Lancés à conditions strictement identiques : comparer un scramblé sans MSA
  à un parent avec MSA ne dirait rien.
- Piège évité : dans l'ordre inverse, les cystéines passent à 23/88/144/218.
  Recopier les contraintes du parent aurait créé des ponts sur les mauvais
  résidus, sans que ça se voie.

## Interprétation
- L'ipSAE est un filtre, pas un classement. On écarte les designs cassés,
  puis on choisit 15 designs divers.
- L'ordre inverse ne diffère que d'un facteur 2,2 en KD (0,82 vs 1,83 nM),
  soit ~0,5 kcal/mol — sous la résolution d'un prédicteur de structure.
  Rapporté, non interprété.
- Le vrai contrôle négatif est le CDR-H3 scramblé.

## Gel
16h30. Dernière séquence à 16h30.

## Amendement après les contrôles — agrégation

Les neuf runs de contrôle (3 contrôles × 3 seeds × 5 samples, réglages
identiques) invalident l'agrégation par médiane seule.

`cdrh3_scrambled` n'est pas décalé, il est **bimodal** : un sample soit
s'arrime vers 0,87, soit s'effondre sous 0,3. Aucun des 45 samples ne tombe
entre 0,42 et 0,75. Les effondrements sont réels, pas un artefact de scoring :
l'ipTM_d0chn suit, n0res passe de 242 à 0, l'interface HER2 de 142 à 0, alors
que le pLDDT moyen reste à 82-85 — les deux chaînes se replient, seule
l'interface échoue.

| contrôle | arrimés | médiane | écart intra-seed |
|---|---|---|---|
| parent_wt | 15/15 | 0,9078 | 0,029 – 0,130 |
| reverse_vl_vh | 15/15 | 0,8842 | 0,061 – 0,121 |
| cdrh3_scrambled | 7/15 | 0,2754 | 0,832 – 0,847 |

Trois conséquences :

- La médiane d'une variable bimodale rapporte le mode qui a gagné les cinq
  tirages. `cdrh3_scrambled_seed1` médiane à 0,8873 et passerait n'importe
  quel filtre de magnitude ; `seed2` à 0,1240. Même séquence, mêmes réglages,
  seed différente. Un design jugé sur une seule seed peut être retenu par
  pur tirage.
- La magnitude d'ipSAE ne sépare pas le négatif : quand le scramblé s'arrime,
  il marque 0,8729 contre 0,9050 pour le parent.
- Ce qui sépare, c'est le **taux d'arrimage** (15/15 contre 7/15) et
  l'**écart intra-seed**, parfaitement disjoints : parent et reverse à
  0,1301 au plus, scramblé à 0,8321 au moins, facteur 6,4.

**Nouvelle règle.** On filtre sur le taux d'arrimage et l'écart, pas sur la
médiane :

1. tous les samples doivent s'arrimer (ipSAE ≥ 0,5) ;
2. écart intra-run ≤ 0,30 ;
3. la médiane n'intervient qu'ensuite, et seulement entre designs qui
   passent — un petit écart entre deux designs retenus ne devient pas pour
   autant une preuve d'affinité.

Seuils calibrés sur les contrôles : les 6 bons runs passent, les 3 scramblés
échouent, chacun sur les deux critères indépendamment. Ils reposent sur n=3
par contrôle et restent provisoires : rien n'établit le taux d'effondrement
réel du parent au-delà de « 0 sur 15 ». À revoir si des designs se groupent
près d'une borne.

L'intuition de garder l'écart comme signal de fiabilité est confirmée. C'est
le meilleur discriminateur des trois ; il devient une barrière, pas un
départage.
