# Searches tree on a network of SGI's running the mpich implementation
# of MPI

CC       =  cc
# Additional CFLAGS: -DSTATS, -DDEBUG
CFLAGS   =  -g -fullwarn -DSTATS
#CFLAGS   =  -g -fullwarn
LDFLAGS  =
INCLUDE  =  -I/usr/local/mpich/include -I.
LIB      =  -L/usr/local/mpich/lib/IRIX/ch_p4 -lmpi

SOURCE = \
        main.c \
	par_tree_search.c \
	par_dfs.c \
	service_requests.c \
	work_remains.c \
	node_stack.c \
	solution.c \
	queue.c \
	terminate.c \
	stats.c \
	cio.c \
	vsscanf.c

HDRS = \
	par_tree_search.h \
	par_dfs.h \
	service_requests.h \
	work_remains.h \
	node_stack.h \
	solution.h \
	queue.h \
	terminate.h \
        stats.h \
	cio.h \
	vsscanf.h

OBJS = \
        main.o \
	par_tree_search.o \
	par_dfs.o \
	service_requests.o \
	work_remains.o \
	node_stack.o \
	solution.o \
	queue.o \
	terminate.o \
        stats.o \
	cio.o \
	vsscanf.o


tree: $(OBJS)
	$(CC) -o tree $(OBJS) $(INCLUDE) $(LIB)

clean:
	rm -f tree *.o core

main.o: cio.h node_stack.h par_tree_search.h solution.h terminate.h \
	stats.h queue.h

par_tree_search.o: cio.h par_tree_search.h node_stack.h par_dfs.h \
        service_requests.h work_remains.h solution.h stats.h

par_dfs.o: par_dfs.h node_stack.h queue.h solution.h stats.h

service_requests.o: service_requests.h node_stack.h queue.h terminate.h \
	stats.h

work_remains.o: work_remains.h node_stack.h terminate.h queue.h \
	service_requests.h stats.h

terminate.o: terminate.h

node_stack.o: node_stack.h

solution.o: solution.h cio.h

queue.o: queue.h node_stack.h

stats.o: stats.h

cio.o: cio.h vsscanf.h

vsscanf.o: vsscanf.h

.c.o:
	$(CC) -c $(CFLAGS) $*.c $(INCLUDE)
