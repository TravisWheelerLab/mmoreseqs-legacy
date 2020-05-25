#include "KSeqWrapper.h"
#include "kseq.h"
#include "FileUtil.h"
#include "Util.h"
#include "Debug.h"
#include <unistd.h>

namespace KSEQFILE {
    KSEQ_INIT(int, read)
}

KSeqFile::KSeqFile(const char* fileName) {
    file = FileUtil::openFileOrDie(fileName, "r", true);
    seq = (void*) KSEQFILE::kseq_init(fileno(file));
}

bool KSeqFile::ReadEntry() {
    KSEQFILE::kseq_t* s = (KSEQFILE::kseq_t*) seq;
    int result = KSEQFILE::kseq_read(s);
    if (result < 0)
        return false;
    entry.headerOffset = s->headerOffset;
    entry.sequenceOffset = s->sequenceOffset;
    entry.multiline = s->multiline;
    entry.name = s->name;
    entry.comment = s->comment;
    entry.sequence = s->seq;
    entry.qual = s->qual;

    return true;
}

KSeqFile::~KSeqFile() {
    kseq_destroy((KSEQFILE::kseq_t*)seq);
    fclose(file);
}


namespace KSEQSTREAM {
    KSEQ_INIT(int, read)
}

KSeqStream::KSeqStream() {
    seq = (void*) KSEQSTREAM::kseq_init(STDIN_FILENO);
}

bool KSeqStream::ReadEntry() {
    KSEQSTREAM::kseq_t* s = (KSEQSTREAM::kseq_t*) seq;
    int result = KSEQSTREAM::kseq_read(s);
    if (result < 0)
        return false;

    entry.name = s->name;
    entry.comment = s->comment;
    entry.sequence = s->seq;
    entry.qual = s->qual;

    return true;
}

KSeqStream::~KSeqStream() {
    kseq_destroy((KSEQSTREAM::kseq_t*)seq);
}

#ifdef HAVE_ZLIB
namespace KSEQGZIP {
    KSEQ_INIT(gzFile, gzread)
}

KSeqGzip::KSeqGzip(const char* fileName) {
    if(FileUtil::fileExists(fileName) == false) {
        errno = ENOENT;
        perror(fileName);
        EXIT(EXIT_FAILURE);
    }

    file = gzopen(fileName, "r");
    if(file == NULL) {
        perror(fileName); EXIT(EXIT_FAILURE);
    }

    seq = (void*) KSEQGZIP::kseq_init(file);
}

bool KSeqGzip::ReadEntry() {
    KSEQGZIP::kseq_t* s = (KSEQGZIP::kseq_t*) seq;
    int result = KSEQGZIP::kseq_read(s);
    if (result < 0)
        return false;

    entry.name = s->name;
    entry.comment = s->comment;
    entry.sequence = s->seq;
    entry.qual = s->qual;
    entry.headerOffset = -1;
    entry.sequenceOffset = -1;
    entry.multiline = s->multiline;

    return true;
}

KSeqGzip::~KSeqGzip() {
    kseq_destroy((KSEQGZIP::kseq_t*)seq);
    gzclose(file);
}
#endif


#ifdef HAVE_BZLIB
namespace KSEQBZIP {

    KSEQ_INIT(BZFILE *, BZ2_bzread)
}

KSeqBzip::KSeqBzip(const char* fileName) {
    if(FileUtil::fileExists(fileName) == false) {
        errno = ENOENT;
        perror(fileName);
        EXIT(EXIT_FAILURE);
    }
    FILE *fp = FileUtil::openFileOrDie(fileName, "r+b", true);
    int bzError;
    file = BZ2_bzReadOpen(&bzError, fp, 0, 0, NULL, 0);
    if(bzError != 0){
        perror(fileName); EXIT(EXIT_FAILURE);
    }
    seq = (void*) KSEQBZIP::kseq_init(file);
}

bool KSeqBzip::ReadEntry() {
    KSEQBZIP::kseq_t* s = (KSEQBZIP::kseq_t*) seq;
    int result = KSEQBZIP::kseq_read(s);
    if (result < 0)
        return false;

    entry.name = s->name;
    entry.comment = s->comment;
    entry.sequence = s->seq;
    entry.qual = s->qual;
    entry.headerOffset = -1;
    entry.sequenceOffset = -1;
    entry.multiline = s->multiline;

    return true;
}

KSeqBzip::~KSeqBzip() {
    kseq_destroy((KSEQBZIP::kseq_t*)seq);
    int bzError;
    BZ2_bzReadClose(&bzError, file);
}
#endif

KSeqWrapper* KSeqFactory(const char* file) {
    KSeqWrapper* kseq = NULL;
    if( strcmp(file, "stdin") == 0 ){
        kseq = new KSeqStream();
        return kseq;
    }

    if(Util::endsWith(".gz", file) == false && Util::endsWith(".bz2", file) == false ) {
        kseq = new KSeqFile(file);
        return kseq;
    }
#ifdef HAVE_ZLIB
    else if(Util::endsWith(".gz", file) == true) {
        kseq = new KSeqGzip(file);
        return kseq;
    }
#else
    else if(Util::endsWith(".gz", file) == true) {
        Debug(Debug::ERROR) << "MMseqs was not compiled with zlib support. Can not read compressed input!\n";
        EXIT(EXIT_FAILURE);
    }
#endif

#ifdef HAVE_BZLIB
    else if(Util::endsWith(".bz2", file) == true) {
        kseq = new KSeqBzip(file);
        return kseq;
    }
#else
    else if(Util::endsWith(".bz2", file) == true) {
        Debug(Debug::ERROR) << "MMseqs was not compiled with bz2lib support. Can not read compressed input!\n";
        EXIT(EXIT_FAILURE);
    }
#endif

    return kseq;
}

