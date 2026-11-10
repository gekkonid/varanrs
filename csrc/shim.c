/*
 * Minimal implementations of the bcftools helpers that the vendored
 * bcftools/filter.c relies on.
 */
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

/* bcftools' error(): print the message and exit, as bcftools does. Called on
 * malformed expressions and on fatal evaluation errors. */
void error(const char *format, ...)
{
    va_list ap;
    fflush(stdout);
    fputs("varanrs: ", stderr);
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    fflush(stderr);
    exit(1);
}
