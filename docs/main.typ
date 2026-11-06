/*
 * Rubber Article Template - Computational Biology TP Report
 * Biodiversity and Structural Evolution of PETase Enzymes in the Global Oceans
 */

#import "@preview/rubber-article:0.5.2": *

// Layout and styling
#show: article.with(
  cols: none,
  eq-chapterwise: false,
  eq-numbering: "(1)",
  header-display: false,
  header-title: "Biodiversity of enzymes in the global oceans",
  lang: "en",
  page-margins: 2.5cm,
  page-paper: "a4",
)

// Frontmatter
#maketitle(
  title: "Biodiversity of enzymes in the global oceans",
  authors: ("Emil Juul Elberg", "Federica Vellano", "Diego Torres Perez"),
  date: datetime.today().display("[day]. [month repr:long] [year]"),
)

#figure(
  image("../fig/world_map_temperature.png", width: 85%),
) <fig:oceans>

#outline(depth: 2)

#v(1fr)

#heading(level: 2, numbering: none, outlined: false)[Code Availability]
#link("https://github.com/Skimaste/PETase.git")[https://github.com/Skimaste/PETase.git]

#pagebreak()

= Introduction
Plastic pollution is one of the most widespread environmental challenges facing marine ecosystems. In 2016, the bacterium _Ideonella sakaiensis_ 201-F6 was isolated and shown to degrade polyethylene terephthalate (PET) using a specialized esterase named PETase @yoshida2016. 

Since then, PETase-like enzymes have been detected across the global ocean, from surface waters down to the deep sea @Alam2020.09.07.285692. Marine microbial communities inhabit widely varying physical regimes characterized by steep gradients in temperature ($-1.6$ to $30.5$ °C), hydrostatic pressure, and nutrients. Homologues of PETase appear diverse and evolve rapidly, although only some may possess true PET-degrading activity, suggesting a broader enzyme family whose distribution and sequence variation reflect environmental adaptation.

Starting from the reference PETase sequence of _I. sakaiensis_ (290 amino acids), we search for homologues in the Tara Oceans microbial gene catalogue, align them to assess sequence diversity, investigate their biogeographical abundance across ocean samples, and map sequence variability onto the 3D structure of the enzyme. Through this workflow, we address three central questions:

+ *How much sequence diversity exists within this enzyme family?*
+ *Is the abundance and distribution of the homologues related to environmental parameters?*
+ *Does sequence diversity map to specific functional or structural regions of the 3D enzyme?*

= Data & Metagenomic Catalogue
Our analysis is based on metagenomic data collected during the Tara Oceans expedition @sunagawa2015, which systematically sampled marine microbial communities across the global ocean. 

The original catalogue of 46 million non-redundant microbial genes was pre-clustered at $E$-value $<= 10^(-3)$ and 70 % reciprocal coverage into 2.4 million clusters. Each cluster is represented by a single sequence and accompanied by an abundance profile across 179 Tara Oceans stations. Temperature in the metadata ranges from $-1.6$ to $30.5$ °C, showing a polar cluster near $0$ °C and a larger tropical/temperate cluster around $25$ °C. Epipelagic samples dominate the dataset: 83 stations are surface ($<= 10$ m), 52 are at the Deep Chlorophyll Maximum (DCM, $10-150$ m), and only 42 are mesopelagic ($> 150$ m, down to $1 000$ m).

The total sequenced read depth varies from 8.1 to 45.6 million reads per sample. To allow unbiased comparisons across stations, the raw abundance of each cluster $c$ at station $s$ was normalized by the total sample sequencing depth:

$ "Normalized Abundance"_(c, s) = ("Read Counts"_(c, s)) / ("Total Library Size"_s) $

= Computational Methods & Experiments

== Looking for homologues
We queried the compressed catalogue (`OM-RGC_v2.rep_seq`) for sequences similar to the PETase of _I. sakaiensis_ (290 aa) using MMseqs2 @steinegger2017. The sensitivity parameter `-s` determines candidate retention before full Smith-Waterman alignment; higher values discover distant homologues at the cost of execution time. We selected `-s 10.0`, well above the default sensitivity of 6.0, identifying 24 hits with an $E$-value below $10^(-3)$.

#figure(
  image("../fig/table.jpeg", width: 85%),
  caption: [MMseqs2 hits obtained with sensitivity `-s 10.0`. Green: coverage `aln_frac` $> 0.6$; blue: target length similar to query (`lengthRatio` between 0.8 and 1.2).],
) <fig:hits>

Catalogue representatives exhibit a wide length distribution (median length of only 143 aa, compared to the 290 aa query), reflecting numerous fragmented open reading frames from metagenomic assemblies. Consequently, an $E$-value cutoff alone is insufficient: a short, highly identical 40 aa match can achieve an $E$-value below $10^(-5)$ while covering only a single helix.

We therefore calculated the aligned fraction over the shorter sequence:

$ "aln_frac" = (e_q - s_q) / min(L_q, L_t) $

where $s_q, e_q$ are alignment coordinates on the query and $L_q, L_t$ are sequence lengths. We retained hits with $"aln_frac" > 0.60$ (covering at least 60 % of the shorter sequence). This rule effectively filters out short local alignments while preserving well-aligned full-domain and partial-domain homologues, yielding 16 reliable homologues for subsequent multiple sequence alignment and structural mapping.

== Aligning the homologues
A multiple sequence alignment (MSA) of the selected homologues was generated using `mmseqs result2msa` and cleaned of trailing null characters (`tr -d '\0'`). The alignment is anchored to the reference PETase query, so each of its 290 positions corresponds to a residue in the reference protein and its 3D crystal structure.

#figure(
  image("../fig/msa2.png", width: 95%),
  caption: [Multiple sequence alignment of the 16 PETase homologues, ordered by sequence identity with the query (most similar at top). Red dashed lines highlight the catalytic residues (G158, S160, D206, H237).],
) <fig:msa>

As illustrated in @fig:msa, alignment coverage varies across sequence positions due to metagenomic fragmentation. The N-terminal signal peptide (residues 1–29) and C-terminal tail (residues 276–290) are poorly covered, whereas the central catalytic core (positions 86–211) is covered by more than half of the homologues. 

Crucially, the canonical esterase/PETase catalytic motif G-x-S-x-G (positions 158–162) is conserved in 8 of the 9 homologues covering that window. The catalytic triad residues (Ser160, Asp206, and His237) are strictly invariant across all spanning homologues, validating that the retrieved marine sequences belong to the bona fide cutinase/PETase family.

== Global Distribution & Oceanographic Niches
We extracted the distribution of PETase homologues across the 179 Tara Oceans stations from the compressed abundance tables. Total PETase abundance was quantified, and the dominant homologue cluster (highest abundance) was determined for each station.

#figure(
  image("../fig/petase_dominant_clusters.png", width: 92%),
  caption: [Global distribution of dominant PETase homologue clusters across 179 Tara Oceans stations on a Robinson projection. Point size is proportional to total PETase abundance at each station.],
) <fig:worldmap>

Rather than a uniform marine distribution, PETase homologues exhibit striking biogeographical partitioning (@fig:worldmap). Four primary clusters account for dominance across 175 of the 179 stations:
- *Cluster `007574371` ($n=61$ stations):* Dominates cold polar and subpolar waters (mean temperature $6.8 plus.minus 8.5$ °C, ranging from $-62.2$°S to $+79.3$°N).
- *Cluster `004492421` ($n=60$ stations):* Dominates warm temperate and tropical epipelagic waters (mean temperature $22.5 plus.minus 5.7$ °C).
- *Cluster `021954915` ($n=30$ stations):* Restricted to warm, shallow surface waters (mean depth 24 m, mean temperature $21.0 plus.minus 4.3$ °C).
- *Cluster `013456435` ($n=24$ stations):* Specialized to the deep mesopelagic zone (mean depth 601 m, mean temperature $9.1 plus.minus 4.8$ °C).

#figure(
  image("../fig/petase_environmental_niches.png", width: 95%),
  caption: [Oceanographic niche differentiation of marine PETase homologues. (A) Thermal and depth niche partitioning of dominant clusters on a log-depth scale. (B) Vertical stratification of total PETase abundance across ocean depth layers (Kruskal-Wallis $H = 22.35, p = 1.40 times 10^(-5)$).],
) <fig:niches>

To assess environmental determinants statistically, we evaluated rank correlations between total PETase abundance and environmental variables:
- Water temperature displays a significant negative correlation (Spearman $rho = -0.344, p = 2.85 times 10^(-6)$), driven by elevated abundances in cold, nutrient-rich polar zones.
- Hydrostatic depth exhibits an inverse correlation (Spearman $rho = -0.183, p = 0.015$).
- Salinity shows a modest negative correlation (Spearman $rho = -0.250, p = 9.08 times 10^(-4)$).

Vertical stratification is pronounced (@fig:niches#text(fill: luma(100))[B]): total abundance peaks at the Deep Chlorophyll Maximum (DCM; mean $914.3 plus.minus 580.4$, median 694.0) and Surface layer (mean $771.2 plus.minus 579.6$, median 595.4), but declines sharply in the Mesopelagic zone (mean $479.3 plus.minus 351.5$, median 342.4). A non-parametric Kruskal-Wallis test confirms that these vertical differences are highly significant ($H = 22.35, p = 1.40 times 10^(-5)$).

= Mapping Diversity on the 3D Structure
To determine whether sequence variation is structurally constrained, we mapped diversity onto the 1.58 Å crystal structure of _I. sakaiensis_ PETase (PDB 6QGC, chain A) @austin2018.

Pairwise alignment between the MSA reference row and the PDB sequence achieved 100 % identity across all 262 resolved residues. For each alignment position $j$, we computed the Shannon information entropy:

$ H_j = - sum_(a in cal(A)) p_j (a) log_2 p_j (a) $

where $p_j(a)$ is the relative frequency of amino acid $a$ at column $j$. For 16 homologues, the maximum attainable entropy is $log_2 16 = 4.00$ bits. To avoid sampling artifacts from fragmented sequences, entropy was evaluated only on columns with $>= 50$ % homologue coverage ($N = 192$ residues). Residue entropies were assigned to the B-factor column of the structure, with unassigned positions set to $-1.0$.

#figure(
  image("../fig/entropy.png", width: 90%),
  caption: [Shannon entropy per position (bars) and sequence coverage $C_j$ (orange curve). Grey bars: positions covered by fewer than 50 % of homologues. Red dashed lines: catalytic triad residues.],
) <fig:entropy>

The entropy profile (@fig:entropy) demonstrates that the catalytic triad (Ser160, Asp206, His237) is strictly invariant ($H_j = 0.00$ bits). In contrast, adjacent positions fluctuate markedly, reaching values above 3.0 bits.

#figure(
  image("../fig/prot.png", width: 50%),
  caption: [Shannon entropy mapped onto the 3D crystal structure of PETase (PDB 6QGC). Blue: conserved; red: variable; grey: unassigned terminal residues. Green spheres highlight the catalytic triad (Ser160, Asp206, His237).],
) <fig:structure>

Visual inspection of the 3D structure (@fig:structure) reveals a clear biophysical pattern: the central $beta$-sheet core of the $alpha$/$beta$-hydrolase fold is strictly conserved (blue), while outer loops and peripheral helices exhibit high evolutionary variability (red).

To quantify this observation, we measured structural packing density as the mean Euclidean distance from each residue's $C_beta$ atom to its 8 nearest spatial $C_beta$ neighbours @krapp2024. $C_beta$ packing distance correlates strongly with Shannon entropy (Spearman $rho = 0.50, p < 0.001$ via 1,000 permutations). In comparison, radial distance from the protein geometric centroid yields a weaker correlation ($rho = 0.33$). This demonstrates that local packing within the folded core—rather than overall spherical depth—is the primary biophysical constraint governing sequence conservation.

== Artifacts and Quality Control
When interpreting metagenomic data, potential computational and experimental artifacts must be carefully examined:

1. *Pre-clustering Bias:* The OM-RGC v2 catalogue clustered sequences at 70 % identity and $E <= 10^(-3)$. While reducing computational complexity, this merges strain-level microdiversity into single representative sequences. The apparent uniformity of each cluster across stations masks internal polymorphism.
2. *Metagenomic Fragmentation:* High environmental complexity leads to fragmented assemblies, reflected by a median gene length of 143 aa in the catalogue. This fragmentation causes drop-offs in MSA coverage at both the N- and C-terminal ends (@fig:msa), restricting reliable entropy calculations to 192 of 290 residues.
3. *Sampling and Depth Bias:* Epipelagic waters ($<= 150$ m) account for 75 % of Tara Oceans stations, and sampling reaches a maximum depth of 1,000 m. True bathypelagic and abyssopelagic communities remain unrepresented. Furthermore, oceanic coverage is uneven, with substantial gaps in the western Pacific.

= Results & Critical Assessment

== Synthesis of Research Questions
Returning to our three motivating questions:
+ *Sequence Diversity:* Marine PETase homologues exhibit substantial sequence divergence (identities from 19.5 % to 46.2 % with _I. sakaiensis_ PETase among filtered homologues, and up to 53.5 % in raw fragments). However, the catalytic G-x-S-x-G motif and active triad are strictly conserved, confirming shared catalytic mechanisms.
+ *Environmental Relationship:* Biogeographic abundance is strongly correlated with temperature ($rho = -0.344$) and depth ($rho = -0.183$). Different homologue clusters occupy distinct thermal and vertical niches, separating into polar surface, warm epipelagic, and deep mesopelagic ecotypes ($p = 1.40 times 10^(-5)$).
+ *3D Structural Localization:* Evolutionary mutations are non-randomly distributed in 3D space: the core $beta$-sheet and catalytic cleft are evolutionarily frozen, whereas solvent-exposed loops accumulate diversity ($rho = 0.50$ with $C_beta$ packing distance).

== Critical Assessment & Limitations
A critical question is whether these marine homologues represent genuine plastic-degrading enzymes in nature:
- *Natural vs. Synthetic Substrates:* PET was introduced into industrial manufacturing only ~70 years ago—an extremely short timescale for microbial speciation. The detected marine enzymes belong to the broader $alpha/beta$-hydrolase and cutinase families, which natively hydrolyze biopolyesters such as algal cutin, suberin, or bacterial polyhydroxyalkanoates (PHA).
- *Homology vs. Biochemical Proof:* While sequence homology and catalytic triad conservation indicate esterase activity, they do not guarantee efficient depolymerization of crystalline PET polymers. Biochemical assays (e.g., clearance zone plates, HPLC analysis of terephthalic acid release) are required before concluding that these enzymes actively degrade oceanic plastic waste.
- *Sample Size Limitations:* After filtering fragmented sequences, our structural conclusions rely on 16 well-covered homologues. While the active site conservation is unequivocal, richer metagenomic databases will be needed to capture subtle co-evolutionary contacts.

=== Contributions
This investigation and report were carried out jointly by Emil Juul Elberg, Federica Vellano, and Diego Torres Perez.

#bibliography("ref.bib")
