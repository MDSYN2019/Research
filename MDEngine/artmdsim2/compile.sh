# script for testing that programs compile correctly (see README)
#!/bin/sh

#set CFLAGS = "-I./src -I/usr/include/openmpi/"
set CFLAGS = "-I./src"

#set LFLAGS = "-L/usr/include/openmpi -lmpi -lpthread -lm"
set LFLAGS = "-lm"

for f in "src/pr_*";
do
  echo $f
  gcc $CFLAGS $f $LFLAGS
done
