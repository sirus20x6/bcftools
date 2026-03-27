// bcf_compat.h -- declarations for C accessor functions for bcf1_t fields
// that contain bitfields or other constructs Zig's cImport cannot handle.
// The implementations are in bcf_compat.c.

#ifndef BCF_COMPAT_H
#define BCF_COMPAT_H

#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>

hts_pos_t bcf_compat_pos(const bcf1_t *rec);
hts_pos_t bcf_compat_rlen(const bcf1_t *rec);
int32_t   bcf_compat_rid(const bcf1_t *rec);
uint16_t  bcf_compat_n_allele(const bcf1_t *rec);
uint16_t  bcf_compat_n_info(const bcf1_t *rec);
uint8_t   bcf_compat_n_fmt(const bcf1_t *rec);
uint32_t  bcf_compat_n_sample(const bcf1_t *rec);
char    **bcf_compat_alleles(const bcf1_t *rec);
char     *bcf_compat_id(const bcf1_t *rec);

// Lifecycle helpers for keeping bcf1_t alive through the pipeline
bcf1_t   *bcf_compat_dup(bcf1_t *rec);
void      bcf_compat_destroy(bcf1_t *rec);

// Header accessor (bcf_hdr_nsamples is a macro)
int       bcf_compat_nsamples(const bcf_hdr_t *hdr);

// Format a bcf1_t record as a VCF text line into a kstring_t
int       bcf_compat_format(bcf_hdr_t *hdr, bcf1_t *rec, kstring_t *str);

// Accessors for bcf_srs_t fields (opaque struct workaround)
bcf_hdr_t *bcf_compat_sr_header(bcf_srs_t *sr, int idx);
int        bcf_compat_sr_next_line(bcf_srs_t *sr);
bcf1_t    *bcf_compat_sr_get_line(bcf_srs_t *sr, int idx);

// Wrappers for htslib macros that cImport may not translate
int bcf_compat_update_info_string(bcf_hdr_t *hdr, bcf1_t *rec, const char *key, const char *val);
int bcf_compat_update_format_int32(bcf_hdr_t *hdr, bcf1_t *rec, const char *key, const int32_t *values, int n);
int bcf_compat_get_genotypes(bcf_hdr_t *hdr, bcf1_t *rec, int32_t **dst, int *ndst);

#endif // BCF_COMPAT_H
