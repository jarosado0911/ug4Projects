mpirun --allow-run-as-root -np 8 ugshell -ex reconstructed_spine_wER.lua -grid wER/Spine5_wER.ugx -numRefs 0 -tstep 5e-6 -endTime 0.015 -outName output -solver ILU
