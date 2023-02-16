#pragma once

#include <stdio.h>
#include <stdlib.h>

static __attribute__((unused)) void pre_info_fn(const char *fname, int line, const char *kind)
{
    fprintf(stderr, "[%s:%d::\033[31;1m%s\033[0m] ", fname, line, kind);
}

#define _usage(msg, kind, terminate)                     \
    do {                                                 \
        pre_info_fn(__FILE__, __LINE__, #kind);          \
        fprintf(stderr, "\033[35;1m%s\033[0m\n", (msg)); \
        if ((terminate)) exit(1);                        \
    } while (0)

#define bug(msg) _usage(msg, bug, 1)
#define err(msg) _usage(msg, error, 1)
#define warn(msg) _usage(msg, warning, 0)
