#pragma once

#include <stddef.h>
#include <sys/types.h>

typedef struct mmapfile *mmapfile_t;

mmapfile_t mmapfile_create_map(char const *fname);
int mmapfile_free_map(mmapfile_t mmf);
ssize_t mmapfile_read(mmapfile_t const mmf, char *buf, size_t pos, size_t count);
size_t mmapfile_filesize(mmapfile_t const mmf);
char const* mmapfile_contents(mmapfile_t const mmf);
