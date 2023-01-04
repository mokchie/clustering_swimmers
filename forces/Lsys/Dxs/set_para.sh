#!/bin/bash

for i in {-4..4}
do
    cd dxs${i}_16
    ./domain.py ${i}
    python set_vars.py 1.0
    python jureca_run.py
    cd ..
done

