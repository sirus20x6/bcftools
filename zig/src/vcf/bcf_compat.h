// bcf_compat.h -- declarations for C accessor functions for bcf1_t fields
// that contain bitfields or other constructs Zig's cImport cannot handle.
// The implementations are in bcf_compat.c.

#ifndef BCF_COMPAT_H
#define BCF_COMPAT_H

#include <htslib/vcf.h>

hts_pos_t bcf_compat_pos(const bcf1_t *rec);
hts_pos_t bcf_compat_rlen(const bcf1_t *rec);
int32_t   bcf_compat_rid(const bcf1_t *rec);
uint16_t  bcf_compat_n_allele(const bcf1_t *rec);
uint16_t  bcf_compat_n_info(const bcf1_t *rec);
uint8_t   bcf_compat_n_fmt(const bcf1_t *rec);
uint32_t  bcf_compat_n_sample(const bcf1_t *rec);
char    **bcf_compat_alleles(const bcf1_t *rec);
char     *bcf_compat_id(const bcf1_t *rec);

#endif // BCF_COMPAT_H
