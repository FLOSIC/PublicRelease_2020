#!/bin/bash

mkdir gnuplot_install
cd gnuplot_install
wget https://downloads.sourceforge.net/project/gnuplot/gnuplot/6.0.1/gnuplot-6.0.1.tar.gz
tar zxvf gnuplot-6.0.1.tar.gz
cd gnuplot-6.0.1/
./configure --prefix=$HOME --without-qt
make
make install

export PATH=$PATH:$HOME/.local/bin:$HOME/bin

echo "PATH=$PATH:$HOME/.local/bin:$HOME/bin" >> $HOME/.bash_profile
echo "export PATH" >> $HOME/.bash_profile


