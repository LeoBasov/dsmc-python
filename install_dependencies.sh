#!/bin/bash

sudo apt update

sudo apt install -y python3
sudo apt install -y python3-pip
sudo apt install -y python3-dbg
    
pip3 install numpy
pip3 install -U pip
pip3 install llvmlite==0.34.0
pip3 install scipy
pip3 install numba
pip3 install pytest
	
echo "INSTALLATION COMPLETE"
