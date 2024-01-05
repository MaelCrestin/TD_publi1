R Notebook
================

Pour des questions pratiques, nous avons pris que 40 séquences (20 seq F
et 20 seq R) importés à partir de la publication : Biodiversity,
environmental drivers, and sustainability of the global deep-sea sponge
microbiome (<https://doi.org/10.1038/s41467-022-32684-4> ; ENA :
PRJNA664762)

Pour charger nos données (nos séquences), on les upload dans nos files,
puis dans le temrinal on utilise “wget -i nom_du_fichier.txt”.

Pour plus de commentaire (voir le tutoriel dada2) :
<https://benjjneb.github.io/dada2/tutorial.html> Chargement du package
dada2.

``` r
library(dada2); packageVersion("dada2")
```

    ## Loading required package: Rcpp

    ## [1] '1.28.0'

Permet de définir la variable de chemin afin qu’elle pointe vers le
répertoire extrait sur la machine.

``` r
path=here::here("Sequence")
list.files(path)
```

    ##  [1] "filtered"                             
    ##  [2] "silva_nr99_v138.1_train_set.fa.gz"    
    ##  [3] "silva_species_assignment_v138.1.fa.gz"
    ##  [4] "SRR12674283_1.fastq.gz"               
    ##  [5] "SRR12674283_2.fastq.gz"               
    ##  [6] "SRR12674290_1.fastq.gz"               
    ##  [7] "SRR12674290_2.fastq.gz"               
    ##  [8] "SRR12674293_1.fastq.gz"               
    ##  [9] "SRR12674293_2.fastq.gz"               
    ## [10] "SRR12674296_1.fastq.gz"               
    ## [11] "SRR12674296_2.fastq.gz"               
    ## [12] "SRR12674298_1.fastq.gz"               
    ## [13] "SRR12674298_2.fastq.gz"               
    ## [14] "SRR12674299_1.fastq.gz"               
    ## [15] "SRR12674299_2.fastq.gz"               
    ## [16] "SRR12674300_1.fastq.gz"               
    ## [17] "SRR12674300_2.fastq.gz"               
    ## [18] "SRR12674302_1.fastq.gz"               
    ## [19] "SRR12674302_2.fastq.gz"               
    ## [20] "SRR12674303_1.fastq.gz"               
    ## [21] "SRR12674303_2.fastq.gz"               
    ## [22] "SRR12674304_1.fastq.gz"               
    ## [23] "SRR12674304_2.fastq.gz"               
    ## [24] "SRR12674306_1.fastq.gz"               
    ## [25] "SRR12674306_2.fastq.gz"               
    ## [26] "SRR12674309_1.fastq.gz"               
    ## [27] "SRR12674309_2.fastq.gz"               
    ## [28] "SRR12674310_1.fastq.gz"               
    ## [29] "SRR12674310_2.fastq.gz"               
    ## [30] "SRR12674313_1.fastq.gz"               
    ## [31] "SRR12674313_2.fastq.gz"               
    ## [32] "SRR12674315_1.fastq.gz"               
    ## [33] "SRR12674315_2.fastq.gz"               
    ## [34] "SRR12674317_1.fastq.gz"               
    ## [35] "SRR12674317_2.fastq.gz"               
    ## [36] "SRR12674318_1.fastq.gz"               
    ## [37] "SRR12674318_2.fastq.gz"               
    ## [38] "SRR12674319_1.fastq.gz"               
    ## [39] "SRR12674319_2.fastq.gz"               
    ## [40] "SRR12674321_1.fastq.gz"               
    ## [41] "SRR12674321_2.fastq.gz"               
    ## [42] "SRR12674322_1.fastq.gz"               
    ## [43] "SRR12674322_2.fastq.gz"

Permet de lire les noms des fichiers “\_1.fastq.gz” et de créer des
liestes correpondantes à des fichiers fastq directs et inverses.

``` r
# Forward and reverse fastq filenames have format "_1.fastq.gz" et "_2.fastq.gz"

fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))

## Extract sample names, assuming filenames have format: .fastq
sample_names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```

Permet de visualiser les profils de qualité des forward read et reverse
read

``` r
plotQualityProfile(fnFs[1:10])
```

![](TD_publi1_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
plotQualityProfile(fnRs[1:10])
```

![](TD_publi1_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

Permet d’attribuer les noms de fichiers fastq.gz filtrés.

``` r
filtFs <- file.path(path, "filtered", paste0(sample_names, "_1.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample_names, "_2.fastq.gz"))
names(filtFs) <- sample_names
names(filtRs) <- sample_names
```

Nous utiliserons les paramètres de filtrage standard : maxN=0(DADA2 ne
nécessite aucun N), truncQ=2, rm.phix=TRUEet maxEE=2. Le maxEEparamètre
définit le nombre maximum d’« erreurs attendues » autorisées dans une
lecture, ce qui constitue un meilleur filtre que la simple moyenne des
scores de qualité.

``` r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)
```

    ##                        reads.in reads.out
    ## SRR12674283_1.fastq.gz    21790     20164
    ## SRR12674290_1.fastq.gz    17462     16319
    ## SRR12674293_1.fastq.gz    26666     25771
    ## SRR12674296_1.fastq.gz    14778     13937
    ## SRR12674298_1.fastq.gz     4545      4357
    ## SRR12674299_1.fastq.gz    36930     35709

Taux d’erreur : L’algorithme DADA2 utilise un modèle d’erreur
paramétrique ( err) et chaque ensemble de données d’amplicon a un
ensemble différent de taux d’erreur. La learnErrorsméthode apprend ce
modèle d’erreur à partir des données, en alternant l’estimation des taux
d’erreur et l’inférence de la composition de l’échantillon jusqu’à ce
qu’ils convergent vers une solution conjointement cohérente. Comme dans
de nombreux problèmes d’apprentissage automatique, l’algorithme doit
commencer par une estimation initiale, pour laquelle les taux d’erreur
maximaux possibles dans ces données sont utilisés (les taux d’erreur si
seule la séquence la plus abondante est correcte et que tous les autres
sont des erreurs).

Erreur Forward :

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 103405680 total bases in 430857 reads from 20 samples will be used for learning the error rates.

Erreur Reverse :

``` r
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 68937120 total bases in 430857 reads from 20 samples will be used for learning the error rates.

Il est toujours intéressant, à titre de contrôle de cohérence, de
visualiser les taux d’erreur estimés :

``` r
plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](TD_publi1_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 20164 reads in 8064 unique sequences.
    ## Sample 2 - 16319 reads in 5038 unique sequences.
    ## Sample 3 - 25771 reads in 5924 unique sequences.
    ## Sample 4 - 13937 reads in 6856 unique sequences.
    ## Sample 5 - 4357 reads in 2365 unique sequences.
    ## Sample 6 - 35709 reads in 8661 unique sequences.
    ## Sample 7 - 19139 reads in 2769 unique sequences.
    ## Sample 8 - 37636 reads in 7395 unique sequences.
    ## Sample 9 - 35566 reads in 7039 unique sequences.
    ## Sample 10 - 17969 reads in 3672 unique sequences.
    ## Sample 11 - 18546 reads in 3814 unique sequences.
    ## Sample 12 - 20760 reads in 8627 unique sequences.
    ## Sample 13 - 14784 reads in 2321 unique sequences.
    ## Sample 14 - 27281 reads in 8588 unique sequences.
    ## Sample 15 - 20039 reads in 4627 unique sequences.
    ## Sample 16 - 25045 reads in 6513 unique sequences.
    ## Sample 17 - 14442 reads in 6087 unique sequences.
    ## Sample 18 - 19155 reads in 3440 unique sequences.
    ## Sample 19 - 19015 reads in 8894 unique sequences.
    ## Sample 20 - 25223 reads in 6478 unique sequences.

``` r
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 20164 reads in 10238 unique sequences.
    ## Sample 2 - 16319 reads in 6727 unique sequences.
    ## Sample 3 - 25771 reads in 7898 unique sequences.
    ## Sample 4 - 13937 reads in 7529 unique sequences.
    ## Sample 5 - 4357 reads in 2843 unique sequences.
    ## Sample 6 - 35709 reads in 10153 unique sequences.
    ## Sample 7 - 19139 reads in 3590 unique sequences.
    ## Sample 8 - 37636 reads in 8987 unique sequences.
    ## Sample 9 - 35566 reads in 8397 unique sequences.
    ## Sample 10 - 17969 reads in 4926 unique sequences.
    ## Sample 11 - 18546 reads in 4860 unique sequences.
    ## Sample 12 - 20760 reads in 10134 unique sequences.
    ## Sample 13 - 14784 reads in 2872 unique sequences.
    ## Sample 14 - 27281 reads in 9332 unique sequences.
    ## Sample 15 - 20039 reads in 5525 unique sequences.
    ## Sample 16 - 25045 reads in 7397 unique sequences.
    ## Sample 17 - 14442 reads in 7730 unique sequences.
    ## Sample 18 - 19155 reads in 4369 unique sequences.
    ## Sample 19 - 19015 reads in 10571 unique sequences.
    ## Sample 20 - 25223 reads in 8507 unique sequences.

``` r
dadaFs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 289 sequence variants were inferred from 8064 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

``` r
mergers <- dada2::mergePairs(
  dadaF = dadaFs,
  derepF = filtFs,
  dadaR = dadaRs,
  derepR = filtRs,
  maxMismatch = 0,
  verbose = TRUE
)
```

    ## 7 paired-reads (in 2 unique pairings) successfully merged out of 18879 (in 5436 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 15151 (in 2706 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 24259 (in 2745 pairings) input.

    ## 5 paired-reads (in 1 unique pairings) successfully merged out of 12998 (in 4429 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 3886 (in 1219 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 33845 (in 5300 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 17959 (in 462 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 36290 (in 5065 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 34148 (in 4610 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 17128 (in 1911 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 17704 (in 2130 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 19426 (in 5829 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 14300 (in 1605 pairings) input.

    ## 2 paired-reads (in 2 unique pairings) successfully merged out of 25801 (in 5052 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 18935 (in 2221 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 23654 (in 3639 pairings) input.

    ## 3 paired-reads (in 1 unique pairings) successfully merged out of 13412 (in 4171 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 18604 (in 3063 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 17761 (in 6748 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 23194 (in 3019 pairings) input.

``` r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

    ## 7 paired-reads (in 2 unique pairings) successfully merged out of 18879 (in 5436 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 15151 (in 2706 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 24259 (in 2745 pairings) input.

    ## 5 paired-reads (in 1 unique pairings) successfully merged out of 12998 (in 4429 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 3886 (in 1219 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 33845 (in 5300 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 17959 (in 462 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 36290 (in 5065 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 34148 (in 4610 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 17128 (in 1911 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 17704 (in 2130 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 19426 (in 5829 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 14300 (in 1605 pairings) input.

    ## 2 paired-reads (in 2 unique pairings) successfully merged out of 25801 (in 5052 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 18935 (in 2221 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 23654 (in 3639 pairings) input.

    ## 3 paired-reads (in 1 unique pairings) successfully merged out of 13412 (in 4171 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 18604 (in 3063 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 17761 (in 6748 pairings) input.

    ## 0 paired-reads (in 0 unique pairings) successfully merged out of 23194 (in 3019 pairings) input.

``` r
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
```

    ##                                                                                                                                                                                                                                                                                                                                                                                                  sequence
    ## 385                           CGACCTACGGGAGGCAGCAGCAGGGAATCTTGCGCAATGGGCGAAAGCCTGACGCAGCGACACCGCGTGGGGGAAGAAGGCCTTAGGGTTGTAAACCCCTTTTCTGTGGGAAGAGATAGGACGGTACCACAGGAATAAGCATCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGATGCAAGCGTTGTCCGGATTTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGAAACCCTGGTAGTCCAGAGAA
    ## 4952 CGACCTACGGGAGGCAGCAGTAAGGAATTGTTCGCAATGGGCGCAAGCCTGACGACGCAACGCCGCGTGGAGGACGAAGATCTTCGGGTTGTAAACTCCTGTCAAGCGGGAAGAACAGCATGCGGGTTAATACCCCGTGTGCCTGACGGTACCGTTGAAGGAAGCCCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGGGGGGGCAAGCGTTGTTCGGAATTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGAAACCCTGGTAGTCCAGAGAA
    ##      abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 385          6       4      78     37         0      0      1   TRUE
    ## 4952         1      34      78     12         0      0      1   TRUE

``` r
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

    ## [1] 20  6

``` r
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
```

    ## 
    ## 331 358 361 363 388 
    ##   2   1   1   1   1

``` r
seqtab
```

    ##             CGACCTACGGGAGGCAGCAGCAGGGAATCTTGCGCAATGGGCGAAAGCCTGACGCAGCGACACCGCGTGGGGGAAGAAGGCCTTAGGGTTGTAAACCCCTTTTCTGTGGGAAGAGATAGGACGGTACCACAGGAATAAGCATCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGATGCAAGCGTTGTCCGGATTTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGAAACCCTGGTAGTCCAGAGAA
    ## SRR12674283                                                                                                                                                                                                                                                                                                                                                                           6
    ## SRR12674290                                                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674293                                                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674296                                                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674298                                                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674299                                                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674300                                                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674302                                                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674303                                                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674304                                                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674306                                                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674309                                                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674310                                                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674313                                                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674315                                                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674317                                                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674318                                                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674319                                                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674321                                                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674322                                                                                                                                                                                                                                                                                                                                                                           0
    ##             GAGTGGCCTACGGGAGGCAGCAGCAGGGAATCTTGCGCAATGGGCGAAAGCCTGACGCAGCGACACCGCGTGGGGGAAGAAGGCCTTAGGGTTGTAAACCCCTTTTCTGTGGGAAGAGATAGGACGGTACCGCAGGAATAAGCATCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGATGCAAGCGTTGTCCGGATTTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGATACCCGGGTAGTCCT
    ## SRR12674283                                                                                                                                                                                                                                                                                                                                                                         0
    ## SRR12674290                                                                                                                                                                                                                                                                                                                                                                         0
    ## SRR12674293                                                                                                                                                                                                                                                                                                                                                                         0
    ## SRR12674296                                                                                                                                                                                                                                                                                                                                                                         5
    ## SRR12674298                                                                                                                                                                                                                                                                                                                                                                         0
    ## SRR12674299                                                                                                                                                                                                                                                                                                                                                                         0
    ## SRR12674300                                                                                                                                                                                                                                                                                                                                                                         0
    ## SRR12674302                                                                                                                                                                                                                                                                                                                                                                         0
    ## SRR12674303                                                                                                                                                                                                                                                                                                                                                                         0
    ## SRR12674304                                                                                                                                                                                                                                                                                                                                                                         0
    ## SRR12674306                                                                                                                                                                                                                                                                                                                                                                         0
    ## SRR12674309                                                                                                                                                                                                                                                                                                                                                                         0
    ## SRR12674310                                                                                                                                                                                                                                                                                                                                                                         0
    ## SRR12674313                                                                                                                                                                                                                                                                                                                                                                         0
    ## SRR12674315                                                                                                                                                                                                                                                                                                                                                                         0
    ## SRR12674317                                                                                                                                                                                                                                                                                                                                                                         0
    ## SRR12674318                                                                                                                                                                                                                                                                                                                                                                         0
    ## SRR12674319                                                                                                                                                                                                                                                                                                                                                                         0
    ## SRR12674321                                                                                                                                                                                                                                                                                                                                                                         0
    ## SRR12674322                                                                                                                                                                                                                                                                                                                                                                         0
    ##             CGACCTACGGGAGGCAGCAGCAGGGAATCTTGCGCAATGGGCGAAAGCCTGACGCAGCGACACCGCGTGGGGGAAGAAGGCCTTAGGGTTGTAAACCCCTTTTCTGTGGGAAGAGATAGGACGGTACCACAGGAATAAGCATCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGATGCAAGCGTTGTCCGGATTTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGATACCCGGGTAGTCCT
    ## SRR12674283                                                                                                                                                                                                                                                                                                                                                                      0
    ## SRR12674290                                                                                                                                                                                                                                                                                                                                                                      0
    ## SRR12674293                                                                                                                                                                                                                                                                                                                                                                      0
    ## SRR12674296                                                                                                                                                                                                                                                                                                                                                                      0
    ## SRR12674298                                                                                                                                                                                                                                                                                                                                                                      0
    ## SRR12674299                                                                                                                                                                                                                                                                                                                                                                      0
    ## SRR12674300                                                                                                                                                                                                                                                                                                                                                                      0
    ## SRR12674302                                                                                                                                                                                                                                                                                                                                                                      0
    ## SRR12674303                                                                                                                                                                                                                                                                                                                                                                      0
    ## SRR12674304                                                                                                                                                                                                                                                                                                                                                                      0
    ## SRR12674306                                                                                                                                                                                                                                                                                                                                                                      0
    ## SRR12674309                                                                                                                                                                                                                                                                                                                                                                      0
    ## SRR12674310                                                                                                                                                                                                                                                                                                                                                                      0
    ## SRR12674313                                                                                                                                                                                                                                                                                                                                                                      0
    ## SRR12674315                                                                                                                                                                                                                                                                                                                                                                      0
    ## SRR12674317                                                                                                                                                                                                                                                                                                                                                                      0
    ## SRR12674318                                                                                                                                                                                                                                                                                                                                                                      3
    ## SRR12674319                                                                                                                                                                                                                                                                                                                                                                      0
    ## SRR12674321                                                                                                                                                                                                                                                                                                                                                                      0
    ## SRR12674322                                                                                                                                                                                                                                                                                                                                                                      0
    ##             CGACCTACGGGAGGCAGCAGTAAGGAATTGTTCGCAATGGGCGCAAGCCTGACGACGCAACGCCGCGTGGAGGACGAAGATCTTCGGGTTGTAAACTCCTGTCAAGCGGGAAGAACAGCATGCGGGTTAATACCCCGTGTGCCTGACGGTACCGTTGAAGGAAGCCCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGGGGGGGCAAGCGTTGTTCGGAATTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGAAACCCTGGTAGTCCAGAGAA
    ## SRR12674283                                                                                                                                                                                                                                                                                                                                                                                                    1
    ## SRR12674290                                                                                                                                                                                                                                                                                                                                                                                                    0
    ## SRR12674293                                                                                                                                                                                                                                                                                                                                                                                                    0
    ## SRR12674296                                                                                                                                                                                                                                                                                                                                                                                                    0
    ## SRR12674298                                                                                                                                                                                                                                                                                                                                                                                                    0
    ## SRR12674299                                                                                                                                                                                                                                                                                                                                                                                                    0
    ## SRR12674300                                                                                                                                                                                                                                                                                                                                                                                                    0
    ## SRR12674302                                                                                                                                                                                                                                                                                                                                                                                                    0
    ## SRR12674303                                                                                                                                                                                                                                                                                                                                                                                                    0
    ## SRR12674304                                                                                                                                                                                                                                                                                                                                                                                                    0
    ## SRR12674306                                                                                                                                                                                                                                                                                                                                                                                                    0
    ## SRR12674309                                                                                                                                                                                                                                                                                                                                                                                                    0
    ## SRR12674310                                                                                                                                                                                                                                                                                                                                                                                                    0
    ## SRR12674313                                                                                                                                                                                                                                                                                                                                                                                                    0
    ## SRR12674315                                                                                                                                                                                                                                                                                                                                                                                                    0
    ## SRR12674317                                                                                                                                                                                                                                                                                                                                                                                                    0
    ## SRR12674318                                                                                                                                                                                                                                                                                                                                                                                                    0
    ## SRR12674319                                                                                                                                                                                                                                                                                                                                                                                                    0
    ## SRR12674321                                                                                                                                                                                                                                                                                                                                                                                                    0
    ## SRR12674322                                                                                                                                                                                                                                                                                                                                                                                                    0
    ##             GAGTGGCCTACGGGAGGCAGCAGACGCGGTTATACGTAAGGTCTAAGTTTTTATTGGAAAAGGCGTAAAGGGTGTGTAGGCAATCATTAAATTTCTCGAATTAATTCGATGTTATGAACTTCCACCAAACAATCAAAAGGGTTTGATTGTTTGGGAGTTATAAATAATAGGATGGTACGGGGTTAAAGGGGTAAATGGAACTTTTATGTAGCGGTAGAATGCTTTAATATGAAATGGAACACTGGAGGCGAAAGCAGTTTACTCGAATCAACCTGACGCTGAAACACGAAAGTAGGGGTAGCAAACAGGATTAGAAACCCTTGTAGTCCGA
    ## SRR12674283                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674290                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674293                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674296                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674298                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674299                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674300                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674302                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674303                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674304                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674306                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674309                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674310                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674313                                                                                                                                                                                                                                                                                                                                           1
    ## SRR12674315                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674317                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674318                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674319                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674321                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674322                                                                                                                                                                                                                                                                                                                                           0
    ##             GAGTGGCCTACGGGAGGCAGCAGACGCGGTTATACGTAAGGTCTAAGTTTTTATTGGAAAAGGCGTAAAGGGTGTGTAGGCAATCATTAAATTTCTCGAATTAATTCGATGTTATGAACTTCCACCAAACAATCAAAAGGGTTTGATTGTTTGGGAGTTATAAATAATAGGATGGTACGGGGTTAAAGGGGTAAATGGAACTTTTATGTAGCGGTAGAATGCTTTAATATGAAATGGAACACTGGAGGCGAAAGCAGTTTACTCGAATCAACCTGACGCTGAAACACGAAAGTAGGGGTAGCAAACAGGATTAGATACCCCGGTAGTCCGA
    ## SRR12674283                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674290                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674293                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674296                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674298                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674299                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674300                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674302                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674303                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674304                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674306                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674309                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674310                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674313                                                                                                                                                                                                                                                                                                                                           1
    ## SRR12674315                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674317                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674318                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674319                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674321                                                                                                                                                                                                                                                                                                                                           0
    ## SRR12674322                                                                                                                                                                                                                                                                                                                                           0

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 0 bimeras out of 6 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1] 20  6

``` r
sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 1

``` r
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample_names
head(track)
```

    ##             input filtered denoisedF denoisedR merged nonchim
    ## SRR12674283 21790    20164     19524     19498      7       7
    ## SRR12674290 17462    16319     15629     15729      0       0
    ## SRR12674293 26666    25771     24933     24929      0       0
    ## SRR12674296 14778    13937     13430     13477      5       5
    ## SRR12674298  4545     4357      4074      4142      0       0
    ## SRR12674299 36930    35709     34558     34819      0       0

``` r
# we save in variable the path to the refdb
# in the working space
silva_train_set <- file.path(path,
                             "silva_nr99_v138.1_train_set.fa.gz")

silva_species_assignment <- file.path(path,
                                      "silva_species_assignment_v138.1.fa.gz")

# then we download the files if they don't already exist

if (!file.exists(silva_train_set)) {
  download.file(
    "https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz",
    silva_train_set,
    quiet = TRUE
  )
}

if (!file.exists(silva_species_assignment)) {
  download.file(
    "https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz",
    silva_species_assignment,
    quiet = TRUE
  )
}
```

``` r
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa
```

    ##                                                                                                                                                                                                                                                                                                                                                                                                      Kingdom   
    ## CGACCTACGGGAGGCAGCAGCAGGGAATCTTGCGCAATGGGCGAAAGCCTGACGCAGCGACACCGCGTGGGGGAAGAAGGCCTTAGGGTTGTAAACCCCTTTTCTGTGGGAAGAGATAGGACGGTACCACAGGAATAAGCATCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGATGCAAGCGTTGTCCGGATTTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGAAACCCTGGTAGTCCAGAGAA                          "Bacteria"
    ## GAGTGGCCTACGGGAGGCAGCAGCAGGGAATCTTGCGCAATGGGCGAAAGCCTGACGCAGCGACACCGCGTGGGGGAAGAAGGCCTTAGGGTTGTAAACCCCTTTTCTGTGGGAAGAGATAGGACGGTACCGCAGGAATAAGCATCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGATGCAAGCGTTGTCCGGATTTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGATACCCGGGTAGTCCT                            "Bacteria"
    ## CGACCTACGGGAGGCAGCAGCAGGGAATCTTGCGCAATGGGCGAAAGCCTGACGCAGCGACACCGCGTGGGGGAAGAAGGCCTTAGGGTTGTAAACCCCTTTTCTGTGGGAAGAGATAGGACGGTACCACAGGAATAAGCATCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGATGCAAGCGTTGTCCGGATTTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGATACCCGGGTAGTCCT                               "Bacteria"
    ## CGACCTACGGGAGGCAGCAGTAAGGAATTGTTCGCAATGGGCGCAAGCCTGACGACGCAACGCCGCGTGGAGGACGAAGATCTTCGGGTTGTAAACTCCTGTCAAGCGGGAAGAACAGCATGCGGGTTAATACCCCGTGTGCCTGACGGTACCGTTGAAGGAAGCCCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGGGGGGGCAAGCGTTGTTCGGAATTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGAAACCCTGGTAGTCCAGAGAA "Bacteria"
    ## GAGTGGCCTACGGGAGGCAGCAGACGCGGTTATACGTAAGGTCTAAGTTTTTATTGGAAAAGGCGTAAAGGGTGTGTAGGCAATCATTAAATTTCTCGAATTAATTCGATGTTATGAACTTCCACCAAACAATCAAAAGGGTTTGATTGTTTGGGAGTTATAAATAATAGGATGGTACGGGGTTAAAGGGGTAAATGGAACTTTTATGTAGCGGTAGAATGCTTTAATATGAAATGGAACACTGGAGGCGAAAGCAGTTTACTCGAATCAACCTGACGCTGAAACACGAAAGTAGGGGTAGCAAACAGGATTAGAAACCCTTGTAGTCCGA                                                          "Bacteria"
    ## GAGTGGCCTACGGGAGGCAGCAGACGCGGTTATACGTAAGGTCTAAGTTTTTATTGGAAAAGGCGTAAAGGGTGTGTAGGCAATCATTAAATTTCTCGAATTAATTCGATGTTATGAACTTCCACCAAACAATCAAAAGGGTTTGATTGTTTGGGAGTTATAAATAATAGGATGGTACGGGGTTAAAGGGGTAAATGGAACTTTTATGTAGCGGTAGAATGCTTTAATATGAAATGGAACACTGGAGGCGAAAGCAGTTTACTCGAATCAACCTGACGCTGAAACACGAAAGTAGGGGTAGCAAACAGGATTAGATACCCCGGTAGTCCGA                                                          "Bacteria"
    ##                                                                                                                                                                                                                                                                                                                                                                                                      Phylum          
    ## CGACCTACGGGAGGCAGCAGCAGGGAATCTTGCGCAATGGGCGAAAGCCTGACGCAGCGACACCGCGTGGGGGAAGAAGGCCTTAGGGTTGTAAACCCCTTTTCTGTGGGAAGAGATAGGACGGTACCACAGGAATAAGCATCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGATGCAAGCGTTGTCCGGATTTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGAAACCCTGGTAGTCCAGAGAA                          "Chloroflexi"   
    ## GAGTGGCCTACGGGAGGCAGCAGCAGGGAATCTTGCGCAATGGGCGAAAGCCTGACGCAGCGACACCGCGTGGGGGAAGAAGGCCTTAGGGTTGTAAACCCCTTTTCTGTGGGAAGAGATAGGACGGTACCGCAGGAATAAGCATCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGATGCAAGCGTTGTCCGGATTTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGATACCCGGGTAGTCCT                            "Chloroflexi"   
    ## CGACCTACGGGAGGCAGCAGCAGGGAATCTTGCGCAATGGGCGAAAGCCTGACGCAGCGACACCGCGTGGGGGAAGAAGGCCTTAGGGTTGTAAACCCCTTTTCTGTGGGAAGAGATAGGACGGTACCACAGGAATAAGCATCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGATGCAAGCGTTGTCCGGATTTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGATACCCGGGTAGTCCT                               "Chloroflexi"   
    ## CGACCTACGGGAGGCAGCAGTAAGGAATTGTTCGCAATGGGCGCAAGCCTGACGACGCAACGCCGCGTGGAGGACGAAGATCTTCGGGTTGTAAACTCCTGTCAAGCGGGAAGAACAGCATGCGGGTTAATACCCCGTGTGCCTGACGGTACCGTTGAAGGAAGCCCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGGGGGGGCAAGCGTTGTTCGGAATTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGAAACCCTGGTAGTCCAGAGAA NA              
    ## GAGTGGCCTACGGGAGGCAGCAGACGCGGTTATACGTAAGGTCTAAGTTTTTATTGGAAAAGGCGTAAAGGGTGTGTAGGCAATCATTAAATTTCTCGAATTAATTCGATGTTATGAACTTCCACCAAACAATCAAAAGGGTTTGATTGTTTGGGAGTTATAAATAATAGGATGGTACGGGGTTAAAGGGGTAAATGGAACTTTTATGTAGCGGTAGAATGCTTTAATATGAAATGGAACACTGGAGGCGAAAGCAGTTTACTCGAATCAACCTGACGCTGAAACACGAAAGTAGGGGTAGCAAACAGGATTAGAAACCCTTGTAGTCCGA                                                          "Proteobacteria"
    ## GAGTGGCCTACGGGAGGCAGCAGACGCGGTTATACGTAAGGTCTAAGTTTTTATTGGAAAAGGCGTAAAGGGTGTGTAGGCAATCATTAAATTTCTCGAATTAATTCGATGTTATGAACTTCCACCAAACAATCAAAAGGGTTTGATTGTTTGGGAGTTATAAATAATAGGATGGTACGGGGTTAAAGGGGTAAATGGAACTTTTATGTAGCGGTAGAATGCTTTAATATGAAATGGAACACTGGAGGCGAAAGCAGTTTACTCGAATCAACCTGACGCTGAAACACGAAAGTAGGGGTAGCAAACAGGATTAGATACCCCGGTAGTCCGA                                                          "Proteobacteria"
    ##                                                                                                                                                                                                                                                                                                                                                                                                      Class                
    ## CGACCTACGGGAGGCAGCAGCAGGGAATCTTGCGCAATGGGCGAAAGCCTGACGCAGCGACACCGCGTGGGGGAAGAAGGCCTTAGGGTTGTAAACCCCTTTTCTGTGGGAAGAGATAGGACGGTACCACAGGAATAAGCATCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGATGCAAGCGTTGTCCGGATTTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGAAACCCTGGTAGTCCAGAGAA                          "Dehalococcoidia"    
    ## GAGTGGCCTACGGGAGGCAGCAGCAGGGAATCTTGCGCAATGGGCGAAAGCCTGACGCAGCGACACCGCGTGGGGGAAGAAGGCCTTAGGGTTGTAAACCCCTTTTCTGTGGGAAGAGATAGGACGGTACCGCAGGAATAAGCATCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGATGCAAGCGTTGTCCGGATTTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGATACCCGGGTAGTCCT                            "Dehalococcoidia"    
    ## CGACCTACGGGAGGCAGCAGCAGGGAATCTTGCGCAATGGGCGAAAGCCTGACGCAGCGACACCGCGTGGGGGAAGAAGGCCTTAGGGTTGTAAACCCCTTTTCTGTGGGAAGAGATAGGACGGTACCACAGGAATAAGCATCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGATGCAAGCGTTGTCCGGATTTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGATACCCGGGTAGTCCT                               "Dehalococcoidia"    
    ## CGACCTACGGGAGGCAGCAGTAAGGAATTGTTCGCAATGGGCGCAAGCCTGACGACGCAACGCCGCGTGGAGGACGAAGATCTTCGGGTTGTAAACTCCTGTCAAGCGGGAAGAACAGCATGCGGGTTAATACCCCGTGTGCCTGACGGTACCGTTGAAGGAAGCCCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGGGGGGGCAAGCGTTGTTCGGAATTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGAAACCCTGGTAGTCCAGAGAA NA                   
    ## GAGTGGCCTACGGGAGGCAGCAGACGCGGTTATACGTAAGGTCTAAGTTTTTATTGGAAAAGGCGTAAAGGGTGTGTAGGCAATCATTAAATTTCTCGAATTAATTCGATGTTATGAACTTCCACCAAACAATCAAAAGGGTTTGATTGTTTGGGAGTTATAAATAATAGGATGGTACGGGGTTAAAGGGGTAAATGGAACTTTTATGTAGCGGTAGAATGCTTTAATATGAAATGGAACACTGGAGGCGAAAGCAGTTTACTCGAATCAACCTGACGCTGAAACACGAAAGTAGGGGTAGCAAACAGGATTAGAAACCCTTGTAGTCCGA                                                          "Alphaproteobacteria"
    ## GAGTGGCCTACGGGAGGCAGCAGACGCGGTTATACGTAAGGTCTAAGTTTTTATTGGAAAAGGCGTAAAGGGTGTGTAGGCAATCATTAAATTTCTCGAATTAATTCGATGTTATGAACTTCCACCAAACAATCAAAAGGGTTTGATTGTTTGGGAGTTATAAATAATAGGATGGTACGGGGTTAAAGGGGTAAATGGAACTTTTATGTAGCGGTAGAATGCTTTAATATGAAATGGAACACTGGAGGCGAAAGCAGTTTACTCGAATCAACCTGACGCTGAAACACGAAAGTAGGGGTAGCAAACAGGATTAGATACCCCGGTAGTCCGA                                                          "Alphaproteobacteria"
    ##                                                                                                                                                                                                                                                                                                                                                                                                      Order          
    ## CGACCTACGGGAGGCAGCAGCAGGGAATCTTGCGCAATGGGCGAAAGCCTGACGCAGCGACACCGCGTGGGGGAAGAAGGCCTTAGGGTTGTAAACCCCTTTTCTGTGGGAAGAGATAGGACGGTACCACAGGAATAAGCATCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGATGCAAGCGTTGTCCGGATTTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGAAACCCTGGTAGTCCAGAGAA                          "SAR202 clade" 
    ## GAGTGGCCTACGGGAGGCAGCAGCAGGGAATCTTGCGCAATGGGCGAAAGCCTGACGCAGCGACACCGCGTGGGGGAAGAAGGCCTTAGGGTTGTAAACCCCTTTTCTGTGGGAAGAGATAGGACGGTACCGCAGGAATAAGCATCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGATGCAAGCGTTGTCCGGATTTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGATACCCGGGTAGTCCT                            "SAR202 clade" 
    ## CGACCTACGGGAGGCAGCAGCAGGGAATCTTGCGCAATGGGCGAAAGCCTGACGCAGCGACACCGCGTGGGGGAAGAAGGCCTTAGGGTTGTAAACCCCTTTTCTGTGGGAAGAGATAGGACGGTACCACAGGAATAAGCATCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGATGCAAGCGTTGTCCGGATTTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGATACCCGGGTAGTCCT                               "SAR202 clade" 
    ## CGACCTACGGGAGGCAGCAGTAAGGAATTGTTCGCAATGGGCGCAAGCCTGACGACGCAACGCCGCGTGGAGGACGAAGATCTTCGGGTTGTAAACTCCTGTCAAGCGGGAAGAACAGCATGCGGGTTAATACCCCGTGTGCCTGACGGTACCGTTGAAGGAAGCCCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGGGGGGGCAAGCGTTGTTCGGAATTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGAAACCCTGGTAGTCCAGAGAA NA             
    ## GAGTGGCCTACGGGAGGCAGCAGACGCGGTTATACGTAAGGTCTAAGTTTTTATTGGAAAAGGCGTAAAGGGTGTGTAGGCAATCATTAAATTTCTCGAATTAATTCGATGTTATGAACTTCCACCAAACAATCAAAAGGGTTTGATTGTTTGGGAGTTATAAATAATAGGATGGTACGGGGTTAAAGGGGTAAATGGAACTTTTATGTAGCGGTAGAATGCTTTAATATGAAATGGAACACTGGAGGCGAAAGCAGTTTACTCGAATCAACCTGACGCTGAAACACGAAAGTAGGGGTAGCAAACAGGATTAGAAACCCTTGTAGTCCGA                                                          "Rickettsiales"
    ## GAGTGGCCTACGGGAGGCAGCAGACGCGGTTATACGTAAGGTCTAAGTTTTTATTGGAAAAGGCGTAAAGGGTGTGTAGGCAATCATTAAATTTCTCGAATTAATTCGATGTTATGAACTTCCACCAAACAATCAAAAGGGTTTGATTGTTTGGGAGTTATAAATAATAGGATGGTACGGGGTTAAAGGGGTAAATGGAACTTTTATGTAGCGGTAGAATGCTTTAATATGAAATGGAACACTGGAGGCGAAAGCAGTTTACTCGAATCAACCTGACGCTGAAACACGAAAGTAGGGGTAGCAAACAGGATTAGATACCCCGGTAGTCCGA                                                          "Rickettsiales"
    ##                                                                                                                                                                                                                                                                                                                                                                                                      Family        
    ## CGACCTACGGGAGGCAGCAGCAGGGAATCTTGCGCAATGGGCGAAAGCCTGACGCAGCGACACCGCGTGGGGGAAGAAGGCCTTAGGGTTGTAAACCCCTTTTCTGTGGGAAGAGATAGGACGGTACCACAGGAATAAGCATCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGATGCAAGCGTTGTCCGGATTTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGAAACCCTGGTAGTCCAGAGAA                          NA            
    ## GAGTGGCCTACGGGAGGCAGCAGCAGGGAATCTTGCGCAATGGGCGAAAGCCTGACGCAGCGACACCGCGTGGGGGAAGAAGGCCTTAGGGTTGTAAACCCCTTTTCTGTGGGAAGAGATAGGACGGTACCGCAGGAATAAGCATCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGATGCAAGCGTTGTCCGGATTTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGATACCCGGGTAGTCCT                            NA            
    ## CGACCTACGGGAGGCAGCAGCAGGGAATCTTGCGCAATGGGCGAAAGCCTGACGCAGCGACACCGCGTGGGGGAAGAAGGCCTTAGGGTTGTAAACCCCTTTTCTGTGGGAAGAGATAGGACGGTACCACAGGAATAAGCATCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGATGCAAGCGTTGTCCGGATTTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGATACCCGGGTAGTCCT                               NA            
    ## CGACCTACGGGAGGCAGCAGTAAGGAATTGTTCGCAATGGGCGCAAGCCTGACGACGCAACGCCGCGTGGAGGACGAAGATCTTCGGGTTGTAAACTCCTGTCAAGCGGGAAGAACAGCATGCGGGTTAATACCCCGTGTGCCTGACGGTACCGTTGAAGGAAGCCCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGGGGGGGCAAGCGTTGTTCGGAATTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGAAACCCTGGTAGTCCAGAGAA NA            
    ## GAGTGGCCTACGGGAGGCAGCAGACGCGGTTATACGTAAGGTCTAAGTTTTTATTGGAAAAGGCGTAAAGGGTGTGTAGGCAATCATTAAATTTCTCGAATTAATTCGATGTTATGAACTTCCACCAAACAATCAAAAGGGTTTGATTGTTTGGGAGTTATAAATAATAGGATGGTACGGGGTTAAAGGGGTAAATGGAACTTTTATGTAGCGGTAGAATGCTTTAATATGAAATGGAACACTGGAGGCGAAAGCAGTTTACTCGAATCAACCTGACGCTGAAACACGAAAGTAGGGGTAGCAAACAGGATTAGAAACCCTTGTAGTCCGA                                                          "Mitochondria"
    ## GAGTGGCCTACGGGAGGCAGCAGACGCGGTTATACGTAAGGTCTAAGTTTTTATTGGAAAAGGCGTAAAGGGTGTGTAGGCAATCATTAAATTTCTCGAATTAATTCGATGTTATGAACTTCCACCAAACAATCAAAAGGGTTTGATTGTTTGGGAGTTATAAATAATAGGATGGTACGGGGTTAAAGGGGTAAATGGAACTTTTATGTAGCGGTAGAATGCTTTAATATGAAATGGAACACTGGAGGCGAAAGCAGTTTACTCGAATCAACCTGACGCTGAAACACGAAAGTAGGGGTAGCAAACAGGATTAGATACCCCGGTAGTCCGA                                                          "Mitochondria"
    ##                                                                                                                                                                                                                                                                                                                                                                                                      Genus
    ## CGACCTACGGGAGGCAGCAGCAGGGAATCTTGCGCAATGGGCGAAAGCCTGACGCAGCGACACCGCGTGGGGGAAGAAGGCCTTAGGGTTGTAAACCCCTTTTCTGTGGGAAGAGATAGGACGGTACCACAGGAATAAGCATCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGATGCAAGCGTTGTCCGGATTTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGAAACCCTGGTAGTCCAGAGAA                          NA   
    ## GAGTGGCCTACGGGAGGCAGCAGCAGGGAATCTTGCGCAATGGGCGAAAGCCTGACGCAGCGACACCGCGTGGGGGAAGAAGGCCTTAGGGTTGTAAACCCCTTTTCTGTGGGAAGAGATAGGACGGTACCGCAGGAATAAGCATCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGATGCAAGCGTTGTCCGGATTTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGATACCCGGGTAGTCCT                            NA   
    ## CGACCTACGGGAGGCAGCAGCAGGGAATCTTGCGCAATGGGCGAAAGCCTGACGCAGCGACACCGCGTGGGGGAAGAAGGCCTTAGGGTTGTAAACCCCTTTTCTGTGGGAAGAGATAGGACGGTACCACAGGAATAAGCATCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGATGCAAGCGTTGTCCGGATTTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGATACCCGGGTAGTCCT                               NA   
    ## CGACCTACGGGAGGCAGCAGTAAGGAATTGTTCGCAATGGGCGCAAGCCTGACGACGCAACGCCGCGTGGAGGACGAAGATCTTCGGGTTGTAAACTCCTGTCAAGCGGGAAGAACAGCATGCGGGTTAATACCCCGTGTGCCTGACGGTACCGTTGAAGGAAGCCCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGGGGGGGCAAGCGTTGTTCGGAATTACTGGGCGTAAAGGGTGTGTAGGCGGCATGTTAAGTCTCATGTGAAATCTCCCGGCTCAACTGGGAGCGGTCATGGGAAACTGGCAAGCTTGAGGGCAGCAGAGGAAAGCGGAATTCCGGGAGTAGTGGTGGAATGCGTAGAAACCCTGGTAGTCCAGAGAA NA   
    ## GAGTGGCCTACGGGAGGCAGCAGACGCGGTTATACGTAAGGTCTAAGTTTTTATTGGAAAAGGCGTAAAGGGTGTGTAGGCAATCATTAAATTTCTCGAATTAATTCGATGTTATGAACTTCCACCAAACAATCAAAAGGGTTTGATTGTTTGGGAGTTATAAATAATAGGATGGTACGGGGTTAAAGGGGTAAATGGAACTTTTATGTAGCGGTAGAATGCTTTAATATGAAATGGAACACTGGAGGCGAAAGCAGTTTACTCGAATCAACCTGACGCTGAAACACGAAAGTAGGGGTAGCAAACAGGATTAGAAACCCTTGTAGTCCGA                                                          NA   
    ## GAGTGGCCTACGGGAGGCAGCAGACGCGGTTATACGTAAGGTCTAAGTTTTTATTGGAAAAGGCGTAAAGGGTGTGTAGGCAATCATTAAATTTCTCGAATTAATTCGATGTTATGAACTTCCACCAAACAATCAAAAGGGTTTGATTGTTTGGGAGTTATAAATAATAGGATGGTACGGGGTTAAAGGGGTAAATGGAACTTTTATGTAGCGGTAGAATGCTTTAATATGAAATGGAACACTGGAGGCGAAAGCAGTTTACTCGAATCAACCTGACGCTGAAACACGAAAGTAGGGGTAGCAAACAGGATTAGATACCCCGGTAGTCCGA                                                          NA

``` r
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```

    ##      Kingdom    Phylum           Class                 Order          
    ## [1,] "Bacteria" "Chloroflexi"    "Dehalococcoidia"     "SAR202 clade" 
    ## [2,] "Bacteria" "Chloroflexi"    "Dehalococcoidia"     "SAR202 clade" 
    ## [3,] "Bacteria" "Chloroflexi"    "Dehalococcoidia"     "SAR202 clade" 
    ## [4,] "Bacteria" NA               NA                    NA             
    ## [5,] "Bacteria" "Proteobacteria" "Alphaproteobacteria" "Rickettsiales"
    ## [6,] "Bacteria" "Proteobacteria" "Alphaproteobacteria" "Rickettsiales"
    ##      Family         Genus
    ## [1,] NA             NA   
    ## [2,] NA             NA   
    ## [3,] NA             NA   
    ## [4,] NA             NA   
    ## [5,] "Mitochondria" NA   
    ## [6,] "Mitochondria" NA

``` r
unqs.mock <- seqtab.nochim["Sequence1"]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Sequence1
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Sequence1 community.\n")
```

    ## DADA2 inferred 0 sample sequences present in the Sequence1 community.

    mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
    match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
    cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")

Erreur dans le type d’argument.

Ici, nous avons pris que 20 séquences pour montrer l’exemple des
démarches à suivre, mais il faudrait beaucoup plus de temps (temps de
chargement) pour une étude complète (avec des milliers de séquences)…
FIN DE DADA2.
