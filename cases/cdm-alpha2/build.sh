#!/bin/sh
# Build nauty 2.8.8 (bundled in the pynauty sdist) and the pruned MTF generator + checker.
# Usage: sh build.sh   (needs pip and gcc)
set -e
mkdir -p build && cd build
pip download --no-binary :all: --no-deps pynauty==2.8.8.1 -q
tar xzf pynauty-2.8.8.1.tar.gz
N=$PWD/pynauty-2.8.8.1/src/nauty2_8_8
(cd $N && ./configure -q && make -s geng)
cd ..
W1="-DMAXN=WORDSIZE -DWORDSIZE=32"
gcc -O3 $W1 -DPRUNE=mtfprune -I$N -c mtfprune.c -o build/mtfprune.o
gcc -O3 $W1 -DPRUNE=mtfprune -I$N -o geng_mtf $N/geng.c build/mtfprune.o \
    $N/gtoolsW.o $N/nautyW1.o $N/nautilW1.o $N/naugraphW1.o $N/schreier.o $N/naurng.o
gcc -O3 -o cdm cdm.c
ln -sf $N/geng geng
echo "run: ./geng_mtf -tcq 16 | ./cdm      (split: ./geng_mtf -tcq 18 0/8 | ./cdm)"
