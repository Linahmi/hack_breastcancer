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

**La marge sur l'écart est mince et à surveiller.** Le scramblé le moins
mauvais est à 0,3481 pour un seuil à 0,30, soit +0,048 ; `seed1` arrime 5/5
et n'échoue que sur l'écart. Le critère d'arrimage ne porte donc plus seul :
l'écart est devenu décisif. Sur la première version du scramblé (99-109,
plus destructrice) la marge était de +0,53. À recalibrer dès que des designs
se groupent près d'une borne, et ne pas traiter 0,30 comme établi.

L'intuition de garder l'écart comme signal de fiabilité est confirmée. C'est
le meilleur discriminateur des trois ; il devient une barrière, pas un
départage.

Note méthodologique sur le plancher de bruit : Aida rapporte l'écart
parent–scramblé contre l'écart inter-seed **du parent** (0,0603), ce qui
donne 1,78× et « résolu ». Nous rapportons aussi le même écart contre celui
du scramblé (0,3357), ce qui donne 0,3×. Les deux sont justes et répondent à
deux questions différentes : « l'écart dépasse-t-il le bruit d'un design
sain ? » (oui) et « une seule seed suffit-elle à trancher ? » (non). Les
deux chiffres bruts sont identiques entre nos deux dérivations.
