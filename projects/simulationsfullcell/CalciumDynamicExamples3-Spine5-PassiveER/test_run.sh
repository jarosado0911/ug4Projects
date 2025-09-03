#!/bin/sh 

mpirun --allow-run-as-root -np 4 ugshell -ex reconstructed_spine_wER.lua -grid Spine5_wER.ugx -numRefs 0 -caInflux 0.0086 -tstep 5e-6 -endTime 0.015 -outName output -solver GS -setting none -minDef 1e-12 -freq 1
