#pragma once

#include <stddef.h>
#include <sys/types.h>

typedef struct fasta_map *fasta_map_t;

fasta_map_t fasta_map_create(char const *fasta_fname);
int fasta_map_free(fasta_map_t fmap);
int fasta_map_query_seq_by_name(fasta_map_t const fmap, char const *name, char *seq);
int fasta_map_query_seq_by_id(fasta_map_t const fmap, size_t id, char *seq);
int fasta_map_query_seq_length_by_name(fasta_map_t const fmap, char const *name, size_t *len);
int fasta_map_query_seq_length_by_id(fasta_map_t const fmap, size_t id, size_t *len);
int fasta_map_query_name_by_id(fasta_map_t const fmap, size_t id, char *name);
int fasta_map_query_id_by_name(fasta_map_t const fmap, char const *name, size_t *id);
int fasta_map_bounds(fasta_map_t const fmap, size_t *max_seq_len, size_t *max_name_len);
size_t fasta_map_num_seqs(fasta_map_t const fmap);
int fasta_map_stats(fasta_map_t const fmap, size_t *numreads, size_t *totbases, size_t *minlen, size_t *maxlen, double *avglen);

ssize_t loadseqs(fasta_map_t const fmap, char **buf, char ***seqs);
