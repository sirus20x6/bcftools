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

bcf1_t   *bcf_compat_dup(bcf1_t *rec)              { return bcf_dup(rec); }
void      bcf_compat_destroy(bcf1_t *rec)           { bcf_destroy(rec); }
int       bcf_compat_nsamples(const bcf_hdr_t *hdr) { return bcf_hdr_nsamples(hdr); }

bcf_hdr_t *bcf_compat_sr_header(bcf_srs_t *sr, int idx) { return sr->readers[idx].header; }
int        bcf_compat_sr_next_line(bcf_srs_t *sr) { return bcf_sr_next_line(sr); }
bcf1_t    *bcf_compat_sr_get_line(bcf_srs_t *sr, int idx) { return bcf_sr_get_line(sr, idx); }

int bcf_compat_format(bcf_hdr_t *hdr, bcf1_t *rec, kstring_t *str) {
    str->l = 0;
    return vcf_format(hdr, rec, str);
}

int bcf_compat_update_info_string(bcf_hdr_t *hdr, bcf1_t *rec, const char *key, const char *val) {
    return bcf_update_info_string(hdr, rec, key, val);
}

int bcf_compat_update_format_int32(bcf_hdr_t *hdr, bcf1_t *rec, const char *key, const int32_t *values, int n) {
    return bcf_update_format_int32(hdr, rec, key, values, n);
}

int bcf_compat_get_genotypes(bcf_hdr_t *hdr, bcf1_t *rec, int32_t **dst, int *ndst) {
    return bcf_get_genotypes(hdr, rec, dst, ndst);
}
