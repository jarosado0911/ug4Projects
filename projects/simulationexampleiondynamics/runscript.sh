#!/bin/sh
ugshell -ex single_neuron_ionDynamics.lua -grid mouse4_standard.ugx -numRefs 0 -dt 0.00001 -endTime 0.05 -vtk -pstep 0.0001 -outName output -nSyn 50 -avgStart 0.000 -avgDur 0.05



