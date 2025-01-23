#!/usr/bin/env bash

git clone https://github.com/genetronhealth/safesim.git
pushd safesim && ./install-dependencies.sh && make && make deploy && popd
gcc -o splitbam.out splitbam.c -I safesim/ext/htslib-1.11-lowdep/ -lhts -L safesim/ext/htslib-1.11-lowdep/
 
