# Table 2 rescored on all copies

vmeasure.sh (run on DRAGEN /staging/tmp/subfam_vmeasure, 2026-11-01): homogeneity, completeness and V-measure (Rosenberg & Hirschberg 2007) on every copy, singletons as their own groups, next to the old coverage/purity columns. Oracle: SubFam and vsearch reproduce Table 2 exactly (groups, singletons, purity). MeShClust2 d7fdb42 rebuilt here with g++ 8.5 does NOT reproduce the published MeShClust2 rows (e.g. B2 0.60: 280 groups vs 636), so its rows here are not comparable with Table 2. cov_ge10 is a fraction.
