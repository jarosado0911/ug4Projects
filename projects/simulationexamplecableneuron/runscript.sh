#!/bin/sh
ugshell -ex single_neuron_somaInjection.lua -grid pyramidal.ugx -dt 1e-5 -endTime 0.05 -vtk -pstep 1e-4 -outName output
