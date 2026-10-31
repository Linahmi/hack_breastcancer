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
