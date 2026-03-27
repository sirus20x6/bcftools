// bcf_compat.c -- C accessor functions for bcf1_t fields that contain
// bitfields which Zig's cImport cannot translate.

#include "bcf_compat.h"

hts_pos_t bcf_compat_pos(const bcf1_t *rec)        { return rec->pos; }
hts_pos_t bcf_compat_rlen(const bcf1_t *rec)       { return rec->rlen; }
int32_t   bcf_compat_rid(const bcf1_t *rec)         { return rec->rid; }
uint16_t  bcf_compat_n_allele(const bcf1_t *rec)    { return rec->n_allele; }
uint16_t  bcf_compat_n_info(const bcf1_t *rec)      { return rec->n_info; }
uint8_t   bcf_compat_n_fmt(const bcf1_t *rec)       { return rec->n_fmt; }
uint32_t  bcf_compat_n_sample(const bcf1_t *rec)    { return rec->n_sample; }
char    **bcf_compat_alleles(const bcf1_t *rec)     { return rec->d.allele; }
char     *bcf_compat_id(const bcf1_t *rec)          { return rec->d.id; }
