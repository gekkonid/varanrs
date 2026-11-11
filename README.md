# varanrs

VARiant ANalysis in RuSt

(pronounced, *Varanus*, like the [genus of monitor lizards](https://en.wikipedia.org/wiki/Varanus_(Varanus)))


# TL;DR 

This toolkit consists of several tools for VCF/BCF manipulation:

* [`varanrs snpsketch`](#snpsketch): generate a summary of a huge VCF/BCF, reporting sample stats, inter-sample distances, and snp-level summaries in a very fast, parallelised manner. Optionally outputs a subset of variants in an R-freindly CSV

* [`varanrs allelefilter`](#allelefilter): separately filter each ALT allele of a VCF, on MAC/MAF. At multiallelic sites, removes alleles individually, returning to an N-1-allelic site (e.g. drops a very rare third allele, yielding a bialelleic site with missing data where the third allele was).

* [`varanrs filter`](#filter): parallel, bcftools-compatible site filtering: `-i/-e` expressions (bcftools syntax), soft filtering, per-allele MAC/MAF filtering, and on-the-fly tag recomputation, all in one pass.

* [`varanrs fill-tags`](#fill-tags): parallel recomputation of AN, AC, AF, MAF, NS, AC_Het/Hom/Hemi, HWE, ExcHet and F_MISSING from genotypes, matching `bcftools +fill-tags`.

* [`varanrs uppercase-alleles`](#uppercase-alleles): forces REF and ALT columns to upper case to work around bugs in various programs.

* And more coming as soon as I need them...


# Usage


## Shared conventions

### Input formats

- Region-parallel subcommands: `.vcf.gz` or `.bcf` with a `.csi` or `.tbi` index. Output is bgzipped VCF, in the same contig order as the input file, always written with a CSI index (`<output>.csi`, as `bcftools ... --write-index`).  Contig lengths come from the header or `--fai`; records beyond a declared (or missing) length are still found via the index, so no records are lost to a wrong header.
- Streaming subcommands: input VCF or BCF from stdin, compressed or not (auto-detected), output VCF text. Single-threaded, designed to be used in some larger filtering pipeline.

Many tools default to streaming, unless you give an input file, in which case they'll use region-parallel iteration and require an index.

All region-parallel subcommands use the following args for input and computation:

| Flag | Description |
|---|---|
| `<INPUT>` | Positional. Variant file, or omit for stdin where supported |
| `--threads <N>` | Worker thread count (default: all cpus) |
| `--contig <NAME>` | Restrict to requested contig(s) (can be given many times) |
| `--fai <PATH>` | FASTA index supplying contig lengths missing from the header (used for window planning and progress) |



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
| `--info-stats info.csv` | Per-INFO-field histograms of site-level values: `field, bin_mean, count`.|


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
field is set with `--info-hist-bins`. **NB: fields are not recalculated**, so
use `bcftools +fill-tags` first if the genotypes have been filtered/modified
since the info column was calculated.

#### Examples

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
| `--fai <PATH>` | *(none)* | FASTA index supplying contig lengths missing from the header; its order is also used to rank sites in `--genotypes` |



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



### `filter`

Filter sites in parallel. Per record, in this order (provided the respective CLI args are given):

1. `-S/--samples-file FILE`: keep only the samples listed in FILE (one per line, output in that order), or with `^FILE` drop them. As `bcftools view -S`, INFO/AC and INFO/AN are recomputed from the kept genotypes, and removing every sample gives sites-only output. Everything below sees only the kept samples.
2. `-I/--include-genotypes` or `-E/--exclude-genotypes EXPR`: per-sample genotype masking, as `bcftools filter -S . -i/-e`. Samples for which the expression is false (`-I`) or true (`-E`) get their GT set to missing (`./.` for diploids, any ploidy kept), and INFO/AC and INFO/AN are decremented if present. E.g. `-E 'FMT/DP<3'` blanks every genotype with DP of 0, 1 or 2. A site-level expression (no FORMAT fields) masks all samples at failing sites. The FILTER column is not touched, and records without a GT field are left unchanged.
3. `--min-ac`/`--min-af`: per-allele filtering exactly as [`allelefilter`](#allelefilter). Sites with no passing ALT allele are dropped.
4. `--fill-tags TAGS`: recompute INFO tags as [`fill-tags`](#fill-tags) does, so the expression sees up-to-date values.
5. `-i/--include` or `-e/--exclude EXPR`: a bcftools filtering expression. Failing sites are dropped, or with `-s/--soft-filter NAME` kept with NAME in FILTER. As with `bcftools filter`, passing sites with an empty FILTER are set to `PASS`, and `-m +`/`-m x` add to existing FILTERs or reset passing sites to PASS.

```
varanrs filter in.bcf -o out.vcf.gz --threads 16 --fill-tags all -e 'QUAL<30 || MAF<0.01 || HWE<1e-6'
varanrs filter in.bcf -o out.vcf.gz -E 'FMT/DP<3 || FMT/GQ<20' --fill-tags all -e 'F_MISSING>0.2'
varanrs filter in.bcf -o out.vcf.gz --min-ac 3 -i 'TYPE="snp" && INFO/DP>10' -s LowQual
bcftools view -Ou in.bcf | varanrs filter -e 'F_MISSING>0.2' | bcftools view -Oz -o out.vcf.gz
```

Expressions are evaluated by bcftools' own expression engine, so the full [bcftools
filtering-expression syntax](https://samtools.github.io/bcftools/bcftools.html#expressions)
is supported, for both the site expression (`-i`/`-e`) and the genotype-masking expression (`-I`/`-E`).
Note that `-S` here subsets samples, as in `bcftools view`; the equivalent of `bcftools filter -S .` is `-I`/`-E`. Perl functions are not available.

#### USAGE

|Input Mode| Support |
|---|---|
|Streaming| YES |
|Region-parallel| YES |

| Flag | Default | Description |
|---|---|---|
| `<INPUT>` | stdin | Indexed VCF/BCF, or omit for stdin |
| `-o, --output <PATH>` | stdout | Output path (required for indexed input) |
| `-S, --samples-file <FILE>` | *(all)* | Subset samples to those listed (`^FILE`: exclude), before all other steps |
| `-I, --include-genotypes <EXPR>` | *(none)* | Set GT to missing for samples where EXPR is not true |
| `-E, --exclude-genotypes <EXPR>` | *(none)* | Set GT to missing for samples where EXPR is true |
| `-i, --include <EXPR>` | *(none)* | Keep sites where EXPR is true |
| `-e, --exclude <EXPR>` | *(none)* | Drop sites where EXPR is true |
| `-s, --soft-filter <NAME>` | *(none)* | Annotate FILTER instead of dropping; `+` picks a name |
| `-m, --mode <+x>` | *(none)* | `+`: add to existing FILTER; `x`: reset passing sites to PASS |
| `--min-ac <N>` / `--min-af <F>` | *(none)* | Per-allele minimum count/frequency |
| `--fill-tags <TAGS>` | *(none)* | Tags to recompute first (see `fill-tags`) |
| `--threads`, `--contig`, `--fai`, `--window-size` | `--window-size 16384` | As for other region-parallel commands. At most 2 × threads windows are held in memory at once. |


### `fill-tags`

Recompute tags from genotypes, producing the same values as `bcftools
+fill-tags` (including its conventions for haploid, half-missing and
polyploid genotypes). Supported tags: `AN`, `AC`, `AF`, `MAF`, `NS`, `AC_Hom`,
`AC_Het`, `AC_Hemi`, `HWE`, `ExcHet`, `F_MISSING` (INFO), `VAF`, `VAF1`
(FORMAT, from `FORMAT/AD`), or `all` (the default, as in bcftools).
Per-population tags (`-S`), `END`, `TYPE` and custom expressions are not
supported.

```
varanrs fill-tags in.bcf -o out.vcf.gz --threads 16
bcftools view -Ou in.bcf | varanrs fill-tags -t AC,AN,AF | bcftools view -Oz -o out.vcf.gz
```

| Flag | Default | Description |
|---|---|---|
| `<INPUT>` | stdin | Indexed VCF/BCF, or omit for stdin |
| `-o, --output <PATH>` | stdout | Output path (required for indexed input) |
| `-t, --tags <LIST>` | `all` | Comma-separated tags to fill |
| `-d, --drop-missing` | off | Do not count half-missing genotypes (`./1`) as hemizygous |
| `--threads`, `--contig`, `--fai`, `--window-size` | `--window-size 16384` | As for other region-parallel commands. At most 2 × threads windows are held in memory at once. |

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
| `--window-size <BP>` | `16384` | Window size in bp (indexed mode) |


# Install/Build

Precompiled binaries are [available on github releases](https://github.com/gekkonid/varanrs/releases).
If these don't work, you'll have to install from source with cargo: Clone this
repo, then run `cargo build --release` and the binary will be in
`target/release/varanrs`. If you don't have cargo, see
[rustup](https://rustup.rs/). Building compiles htslib from source (statically
linked), so a C compiler and libclang (for bindgen; e.g. `apt install
libclang-dev`) are needed at build time. The resulting binary should only
depend on the C runtime (so you may need to compile from source if you are
using an old Linux release, e.g. on clusters).

