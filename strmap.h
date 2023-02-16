#ifndef STRMAP_H_
#define STRMAP_H_

#include <stddef.h>
#include "uthash.h"

struct keyval
{
    char const *key;
    size_t value;
    UT_hash_handle hh;
};

typedef struct { struct keyval *ht; } strmap_t;

void strmap_init(strmap_t *strmap)
{
    strmap->ht = NULL;
}

void strmap_free(strmap_t *strmap)
{
    struct keyval *k, *tmp;

    HASH_ITER(hh, strmap->ht, k, tmp)
    {
        HASH_DEL(strmap->ht, k);
        free(k);
    }
}

int strmap_insert(strmap_t *strmap, char const *key, size_t keylen, size_t value)
{
    struct keyval *k = malloc(sizeof(*k));
    k->key = key;
    k->value = value;
    HASH_ADD_KEYPTR(hh, strmap->ht, k->key, keylen, k);
    return 0;
}

int strmap_insert2(strmap_t *strmap, char const *key, size_t value)
{
    return strmap_insert(strmap, key, strlen(key), value);
}

int strmap_query(strmap_t *strmap, char const *key, size_t *value)
{
    struct keyval *k;
    HASH_FIND_STR(strmap->ht, key, k);

    if (!k) return -1;

    *value = k->value;
    return 0;
}

#endif
