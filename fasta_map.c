#include "fasta_map.h"
#include "mmapfile.h"
#include "strmap.h"
#include "usage.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <limits.h>
#include <stddef.h>
#include <sys/types.h>


typedef struct
{
    char *nameptr; /* read name (also key for UT_hash) */
    size_t len;    /* read length */
    size_t pos;    /* file offset position */
    size_t ll;     /* FASTA line width for this read */
} faidx_record_t;

struct fasta_map
{
    strmap_t idtab;
    faidx_record_t *records;
    mmapfile_t faidx, fasta;
    size_t numseqs, namebnd, seqbnd;
};

static inline void parse_mapped_faidx(fasta_map_t fmap)
{
    faidx_record_t *recs, *rec;
    size_t n, m;
    n = m = 0;
    recs = NULL;

    ssize_t remain = mmapfile_filesize(fmap->faidx);
    char *ptr = (char*)mmapfile_contents(fmap->faidx);
    char line[4096];

    while (remain > 0)
    {
        char *nl = memchr(ptr, '\n', remain);
        if (!nl) break;
        size_t len = nl - ptr;
        memcpy(line, ptr, len);
        line[len] = (char)0;

        if (n+1 > m)
        {
            m += 256;
            recs = realloc(recs, m * sizeof(faidx_record_t));

            if (!recs)
            {
                err("ran out of memory");
            }
        }

        rec = &recs[n];

        if (sscanf(line, "%*s %lu %lu %lu %*u", &rec->len, &rec->pos, &rec->ll) <= 0)
            err("sscanf faidx");

        rec->nameptr = ptr;
        remain -= (len+1);
        ptr = nl+1;
        size_t namelen = 0;

        while (!isspace(rec->nameptr[namelen]))
            ++namelen;

        strmap_insert(&fmap->idtab, rec->nameptr, namelen, n++);

        fmap->seqbnd = fmap->seqbnd > rec->len? fmap->seqbnd : rec->len;
        fmap->namebnd = fmap->namebnd > namelen? fmap->namebnd : namelen;
    }

    recs = realloc(recs, n * sizeof(faidx_record_t));

    if (!recs) err("realloc release");

    fmap->numseqs = n;
    fmap->records = recs;
}

fasta_map_t fasta_map_create_internal(char const *fasta_fname, char const *faidx_fname)
{
    fasta_map_t fmap = malloc(sizeof(*fmap));
    fmap->faidx = mmapfile_create_map(faidx_fname);
    fmap->fasta = mmapfile_create_map(fasta_fname);
    strmap_init(&fmap->idtab);
    fmap->namebnd = fmap->seqbnd = 0;
    parse_mapped_faidx(fmap);
    return fmap;
}

int fasta_map_free(fasta_map_t fmap)
{
    mmapfile_free_map(fmap->fasta);
    mmapfile_free_map(fmap->faidx);
    strmap_free(&fmap->idtab);
    free(fmap->records);
    memset(fmap, 0, sizeof(*fmap));
    free(fmap);
    return 0;
}

static int get_seq_internal(char const *filecontents, faidx_record_t const *record, char *seq)
{
    size_t ll = record->ll;
    size_t locpos = 0;
    ssize_t remain = record->len;
    char *bufptr = (char*)filecontents + record->pos;
    char *dstptr = seq;

    while (remain > 0)
    {
        size_t cnt = ll < remain? ll : remain;
        memcpy(dstptr, bufptr + locpos, cnt);
        dstptr += cnt;
        remain -= cnt;
        locpos += (cnt+1);
    }

    seq[record->len] = (char)0;
    return 0;
}

int fasta_map_query_seq_by_name(fasta_map_t const fmap, char const *name, char *seq)
{
    size_t id;

    if (fasta_map_query_id_by_name(fmap, name, &id))
        return -1;

    return fasta_map_query_seq_by_id(fmap, id, seq);
}

int fasta_map_query_seq_by_id(fasta_map_t const fmap, size_t id, char *seq)
{
    faidx_record_t record = fmap->records[id];
    return get_seq_internal(mmapfile_contents(fmap->fasta), &record, seq);
}


int fasta_map_query_seq_length_by_name(fasta_map_t const fmap, char const *name, size_t *len)
{
    size_t id;

    if (fasta_map_query_id_by_name(fmap, name, &id))
        return -1;

    return fasta_map_query_seq_length_by_id(fmap, id, len);
}

int fasta_map_query_seq_length_by_id(fasta_map_t const fmap, size_t id, size_t *len)
{
    if (!fmap || !len) return -1;
    *len = fmap->records[id].len;
    return 0;
}

int fasta_map_query_name_by_id(fasta_map_t const fmap, size_t id, char *name)
{
    size_t namelen = 0;
    char *nameptr = fmap->records[id].nameptr;

    while (!isspace(nameptr[namelen]))
        ++namelen;

    memcpy(name, nameptr, namelen);
    name[namelen] = (char)0;
    return 0;
}

int fasta_map_query_id_by_name(fasta_map_t const fmap, char const *name, size_t *id)
{
    size_t idx;

    if (!fmap || !id) return -1;

    if (strmap_query(&fmap->idtab, name, &idx))
        return -1;

    *id = idx;
    return 0;
}

int fasta_map_bounds(fasta_map_t const fmap, size_t *max_seq_len, size_t *max_name_len)
{
    if (!fmap) return -1;
    if (max_seq_len) *max_seq_len = fmap->seqbnd;
    if (max_name_len) *max_name_len = fmap->namebnd;
    return 0;
}

size_t fasta_map_num_seqs(fasta_map_t const fmap)
{
    return fmap->numseqs;
}

int fasta_map_stats(fasta_map_t const fmap, size_t *numreads, size_t *totbases, size_t *minlen, size_t *maxlen, double *avglen)
{
    size_t _numreads = fasta_map_num_seqs(fmap);
    size_t _totbases = 0;
    size_t _maxlen = 0;
    size_t _minlen = SIZE_MAX;

    for (size_t i = 0; i < _numreads; ++i)
    {
        size_t len = (fmap->records + i)->len;
        _totbases += len;
        _maxlen = len > _maxlen? len : _maxlen;
        _minlen = len < _minlen? len : _minlen;
    }

    double _avglen = ((double)_totbases) / ((double)_numreads);

    if (numreads) *numreads = _numreads;
    if (totbases) *totbases = _totbases;
    if (minlen) *minlen = _minlen;
    if (maxlen) *maxlen = _maxlen;
    if (avglen) *avglen = _avglen;

    return 0;
}

fasta_map_t fasta_map_create(char const *fasta_fname)
{
    int len = strlen(fasta_fname);
    char *faidx_fname = alloca(len+5);
    strcpy(faidx_fname, fasta_fname);
    strcat(faidx_fname, ".fai");

    return fasta_map_create_internal(fasta_fname, faidx_fname);
}

ssize_t loadseqs(fasta_map_t const fmap, char **buf, char ***seqs)
{
    if (!fmap || !buf || !seqs)
        return -1;

    size_t buflen = 0;
    size_t numseqs = fasta_map_num_seqs(fmap);

    /*
     * Pre-compute the total number of bytes needed
     * to store sequences and their null-terminals.
     */
    for (size_t i = 0; i < numseqs; ++i)
    {
        size_t len;
        fasta_map_query_seq_length_by_id(fmap, i, &len);
        buflen += (len+1);
    }

    char **_seqs = malloc(numseqs*sizeof(char*));
    char *_buf = malloc(buflen);
    char *ptr = _buf;

    for (size_t i = 0; i < numseqs; ++i)
    {
        size_t len;
        fasta_map_query_seq_length_by_id(fmap, i, &len);
        fasta_map_query_seq_by_id(fmap, i, ptr);
        _seqs[i] = ptr;
        ptr += (len+1);
    }

    *buf = _buf;
    *seqs = _seqs;

    return numseqs;
}
