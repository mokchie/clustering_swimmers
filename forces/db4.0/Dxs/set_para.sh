#!/bin/bash

for i in {0..8}
do
    cd dxs${i}_16
    python set_vars.py 4.0
    python jureca_run.py
    cd ..
done

