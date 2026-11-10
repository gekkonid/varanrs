# Vendored bcftools sources

`filter.c`, `filter.h` and `bcftools.h` are copied from bcftools 1.23
(https://github.com/samtools/bcftools, tag `1.23`), MIT/Expat licensed (see
`LICENSE`). They provide the `-i/--include` / `-e/--exclude` expression engine
used by `varanrs filter`, so expressions behave exactly as in bcftools.

Local changes (keep this list up to date when re-vendoring):

1. `filter.c`: removed a leftover debug `fprintf(stderr,"pass=%d [%s]\n",...)`
   in `cmp_vector_strings()` that bcftools 1.23 prints on every site-level
   regex comparison.

Helpers that `filter.c` expects from the rest of bcftools (`error()`) are
provided by `../shim.c`. Perl filters are disabled (`ENABLE_PERL_FILTERS 0`).
