#include "mmapfile.h"
#include "usage.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>   /* SSIZE_MAX */
#include <fcntl.h>    /* open */
#include <sys/stat.h> /* struct stat */
#include <sys/mman.h> /* mmap */
#include <unistd.h>   /* close */

struct mmapfile
{
    char *filecontents;
    size_t filesize;
};

mmapfile_t mmapfile_create_map(char const *fname)
{
    mmapfile_t mmf = NULL;
    struct stat st;
    size_t filesize;
    void *ptr;
    int fd;

    if (stat(fname, &st))
        err("stat()");

    filesize = st.st_size;

    if ((fd = open(fname, O_RDONLY)) < 0)
        err("open()");

    if (!(mmf = malloc(sizeof(*mmf))))
        goto error;

    if ((ptr = mmap(NULL, filesize, PROT_READ, MAP_PRIVATE, fd, 0)) == MAP_FAILED)
        goto error;

    close(fd);
    mmf->filecontents = (char *)ptr;
    mmf->filesize = filesize;
    return mmf;

error:
    close(fd);
    free(mmf);
    return NULL;
}

int mmapfile_free_map(mmapfile_t mmf)
{
    if (!mmf) return -1;

    void *ptr = (void*)mmf->filecontents;
    munmap(ptr, mmf->filesize);
    memset(mmf, 0, sizeof(*mmf));
    free(mmf);

    return 0;
}

ssize_t mmapfile_read(mmapfile_t const mmf, char *buf, size_t pos, size_t count)
{
    if (!mmf || !buf || pos >= mmf->filesize)
        return -1;

    if (pos + count > mmf->filesize)
    {
        count = mmf->filesize - pos;
        count = count < SSIZE_MAX? count : SSIZE_MAX;
    }

    memcpy(buf, mmf->filecontents + pos, count);
    buf[count] = (char)0;

    return count;
}

size_t mmapfile_filesize(mmapfile_t const mmf)
{
    return mmf? mmf->filesize : 0;
}

char const* mmapfile_contents(mmapfile_t const mmf)
{
    return mmf? mmf->filecontents : NULL;
}
