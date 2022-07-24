#!/bin/bash

echo "Run Docker Test"
# source /root/OpenFOAM/OpenFOAM-v1906/etc/bashrc

echo "Copying input data from tests"
cp -r /test/taylor_green_vortex/* .
cp -r /test/execute-tests.sh .

echo "Running execute script"
bash ./execute-tests.sh || exit 1

echo "Copying to outputs"
# cp -r * /outputs/

# echo 'simulation_status: succeeded' 2>&1 | tee "ucns3d.log"