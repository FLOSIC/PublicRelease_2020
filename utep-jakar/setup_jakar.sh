#!/bin/bash

mkdir lapack_install
cd lapack_install
wget https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.10.0.tar.gz
tar zxvf v3.10.0.tar.gz
cd lapack-3.10.0/
cp make.inc.example make.inc
make

mkdir ~/lib
cp liblapack.a ~/lib/
