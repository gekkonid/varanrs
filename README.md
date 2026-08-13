# varanrs

VARiant ANalysis in RuSt

(pronounced, *Varanus*, like the genus of monitor lizards)


# TL;DR 

This toolkit consists of several tools for VCF/BCF manipulation:

* `varanrs snpsketch`: generate a very fast, parallelised summary of a huge VCF/BCF, reporting sample stats (missingness, het rate, depth), inter-sample distances, and snp-level summaries (histograms over any numeric INFO/XXXX fields). Optionally outputs a subset of variants in an R-freindly CSV
* `varanrs allelefilter`: separately filter each ALT allele of a VCF, on MAC/MAF. At multiallelic sites, removes alleles individually, returning to an N-1-allelic site (e.g. drops a very rare third allele, yielding a bialelleic site with missing data where the third allele was).
* `varanrs uppercase-alleles`: forces REF and ALT columns to upper case to work around bugs in various programs.
* And more coming as soon as I need them...


# Usage


## Shared conventions

### Input formats

- **Region-parallel subcommands**: VCF.gz+`.csi`/`.tbi` or BCF+`.csi`. Index must be exist. Either VCF header MUST have contig with length, or you must pass a fasta.fai with --fai. 
- **Streaming subcommands**: uncompressed VCF text or uncompressed BCF binary (auto-detected). Single-threaded, designed to be used in some larger filtering pipeline.

Many tools default to streaming, unless you give an input file, in which case they'll use region-parallel iteration and require an index.


All region-parallel subcommands use the following args for input and computation:

| Flag | Description |
|---|---|
| `<INPUT>` | Positional. Variant file, or omit for stdin on supported commands |
| `--threads <N>` | Worker thread count (default: all cpus) |
| `--contig <NAME>` | Restrict to contig(s) (can be given many times) |
| `--fai <PATH>` | FASTA index for contig names, lengths, and custom ordering, must be used when the VCF header lacks `##contig` length info |



## Commands

### `snpsketch`

Subsample an indexed VCF/BCF to estimate population-genetic summary statistics.
Reads 1 in every N genomic windows (e.g every 10th 5kb window: 1-5000,
50001-55000, ...) and reports:


| Output file | Contents |
|---|---|
| `--sample-stats sample_stats.csv` | Per-sample stats: `sample_id, n_missing, n_total, miss_rate, n_het, het_rate, avg_dp` |
| `--pairs pairs.csv` | Pairwise inter-sample distances (upper triangle): `sample_i, sample_j, sample_i_id, sample_j_id, n_diff, n_common, distance` |
| `--genotypes out.csv` | Per-site genotype matrix: `0`=AA, `1`=AB, `2`=BB, `.`=missing |
| `--info-stats out.csv` | Per-INFO-field histograms of site-level values: `field, bin_mean, count` |


Distance in `--pairs` is the Hamming distance on called genotypes: `n_diff / (2
* n_obs)`. One can read this in R and do ordination or clustering, e.g. using
`cmdscale()` (see below for snippet).

In `--sample-stats`, `miss_rate` is the fraction of sites with missing
genotypes, `avg_dp` is the mean FORMAT/DP over non-missing sites (empty when no
DP values are present), and `het_rate` is the fraction of called sites that are
heterozygous.

`--info-stats` builds a histogram over sites for every numeric INFO field. It
reports these as a tidy/long format csv, with columns field, bin_mean (mean of
the histogram bin), and count (number of sites). For fields with Number=A or
Number=R, only the first ALT allele's value is used and non-numeric fields or
per-genotype fields (Number=G) are skipped. The maximum number of bins per
field is set with `--info-hist-bins`.

```
varanrs snpsketch huge.bcf --stride 100 --threads 8
varanrs snpsketch huge.bcf --fai ref.fai --genotypes gts.csv
```

#### USAGE

|Input Mode| Support |
|---|---|
|Streaming| NO
|Region-parallel| YES |

| Flag | Default | Description |
|---|---|---|
| `<INPUT>` | *(required)* | Indexed VCF/BCF (`.tbi` or `.csi` alongside) |
| `--stride <N>` | `100` | Process 1 in every N windows. 1 means all windows. |
| `--chunk <BP>` | `16384` | Genomic window size in bp |
| `--pairs <PATH>` | `pairs.csv` | Pairwise distance output |
| `--sample-stats <PATH>` | `sample_stats.csv` | Per-sample missingness, heterozygosity, and depth output |
| `--genotypes <PATH>` | *(none)* | Optional site×sample genotype matrix |
| `--info-stats <PATH>` | *(none)* | Optional per-INFO-field streaming histograms |
| `--info-hist-bins <N>` | `100` | Maximum bins per INFO field histogram |
| `--threads <N>` | system parallelism | Worker threads |
| `--contig <NAME>` | all | Restrict to named contig(s) (repeatable) |
| `--fai <PATH>` | *(none)* | FASTA index for contig names, lengths, and ordering. Can be used to subset chroms. |



To see a PCA (well, MDS), in R:

```
library(tidyverse)
dat = read_csv("pairs.csv") |>
    mutate(sample_i_id=as.character(sample_i_id),sample_j_id=as.character(sample_j_id)) |>
    arrange(sample_i, sample_j)
allsamp = unique(c(dat$sample_i_id, dat$sample_j_id))
d = structure(dat$distance, Size = length(allsamp), Labels = allsamp, Diag = FALSE, Upper = FALSE, class = "dist")
plot(cmdscale(d))
```



### `allelefilter`

Independently filter each allele at a site, by minimum allele count (AC) and/or
allele frequency (AF). This is primarily useful to purge "pseudo-multiallelic"
sites, where there are a tiny number of a third allele at an otherwise valid
biallelic site. 

```
varanrs allelefilter in.vcf -o out.vcf --min-ac 5 --min-af 0.01
bcftools view -Ou huge.bcf | varanrs allelefilter --min-ac 3 | bcftools +fill-tags -Ou -- -t all | bcftools view -Oz -o filtered.vcf.gz
```

#### USAGE

|Input Mode| Support |
|---|---|
|Streaming| YES |
|Region-parallel| YES |

| Flag | Default | Description |
|---|---|---|
| `[INPUT]` | stdin | Input VCF/BCF path, or `-` for stdin |
| `-o, --output <PATH>` | stdout | Output VCF path, or `-` for stdout |
| `--min-ac <N>` | *(none)* | Minimum allele count (inclusive) |
| `--min-af <F>` | *(none)* | Minimum allele frequency (inclusive) |



### `uppercase-alleles`

Force REF and ALT alleles to ASCII-uppercase. Workaround for GLnexus, which
occasionally emits lowercase allele characters that downstream tools reject.


```
varanrs uppercase-alleles in.bcf --output out.bcf.vcf.gz --threads 4
bcftools view -Ou in.bcf | varanrs uppercase-alleles | bcftools view -Oz -o out.vcf.gz
```

#### USAGE

|Input Mode| Support |
|---|---|
|Streaming| YES |
|Region-parallel| YES |



| Flag | Default | Description |
|---|---|---|
| `<INPUT>` | stdin | Indexed VCF/BCF file, or omit for stdin |
| `-o, --output <PATH>` | stdout | Output VCF path (required for indexed input) |
| `--threads <N>` | system parallelism | Worker threads (indexed mode only) |
| `--contig <NAME>` | all | Restrict to named contig(s) |
| `--fai <PATH>` | (none) | FASTA index for contig names/lengths |
| `--window-size <BP>` | `1000000` | Window size in bp (indexed mode) |


# Install/Build


Clone, then 

```
cargo build --release
```

and the binary will be in `target/release/varanrs`
