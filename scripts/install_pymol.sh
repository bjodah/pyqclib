#!/bin/bash
wget http://switch.dl.sourceforge.net/project/pymol/pymol/1.7/pymol-v1.7.2.1.tar.bz2
tar xvjf pymol-v1.7.2.1.tar.bz2
cd pymol*
python setup.py build
python setup.py install --user
