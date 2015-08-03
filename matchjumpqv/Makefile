#CC=gcc
CC=icc
CXX=icpc
CFLAGS=-O3 -ip -ipo 
CXXFLAGS=-O3 -ip -ipo  -D__STDC_FORMAT_MACROS -ltbb 
#CFLAGS= -D_GNU_SOURCE -g -O -Wall
#-D_GNU_SOURCE suppress warning for lround
all: matchjumpqv matchdump match2png readtag
clean:
	rm matchjumpqv readtag matchdump match2png
matchjumpqv: matchjumpqv.c basecache.h colorbase.h seqmap.c seqmap.h colorbase.c colorbase.h  basecache.h basecache.c adapter.c adapter.h
	$(CC) $(CFLAGS)  matchjumpqv.c seqmap.c colorbase.c basecache.c adapter.c -o matchjumpqv -lm
matchdump: matchdump.c basecache.c seqmap.c colorbase.c matchlocations.c matchlocfile.h basecache.h colorbase.h seqmap.h 
	${CC} ${CFLAGS} matchdump.c basecache.c colorbase.c seqmap.c matchlocations.c -o matchdump
match2png: match2png.c basecache.c seqmap.c colorbase.c matchlist.c basecache.h matchlocfile.h colorbase.h seqmap.h matchlist.h view.h matchlocations.h
	${CC} ${CFLAGS} match2png.c basecache.c colorbase.c seqmap.c matchlist.c matchlocations.c -o match2png -lm -lgd
readtag: readtag.c matchlocations.c matchlocfile.c basecache.h colorbase.h matchlocations.h matchlocfile.h seqmap.c
	${CC} ${CFLAGS} readtag.c matchlocations.c matchlocfile.c seqmap.c -o readtag

