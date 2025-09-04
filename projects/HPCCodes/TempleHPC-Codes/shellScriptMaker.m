function shellScriptMaker(nFiles)
nn=linspace(0,nFiles,nFiles+1);
density = -0.01;

for k=1:length(nn)
    j=k-1;
    filename=sprintf('Runs/%s_%i.sh','script',j);
    fileID = fopen(filename,'wt');
    density = density + 0.01;
    fprintf(fileID,'#!/bin/sh \n');
    fprintf(fileID,'#PBS -l walltime=00:60:00 \n');
    fprintf(fileID,sprintf('#PBS -N Cell6RyR-job%i\n',j));
    fprintf(fileID,'#PBS -q normal \n');
    fprintf(fileID,'#PBS -l nodes=1:ppn=28 \n');
    fprintf(fileID,'#PBS -e error.txt \n');
    fprintf(fileID,'#PBS -o output.txt \n');
    
    fprintf(fileID,'cd $PBS_O_WORKDIR \n\n');
    
    %fprintf(fileID,'mpirun -np 28 ugshell -ex /home/tuh36651/work/Spine-Dendrite-Calcium/reconstructed_spine_wER.lua -grid /home/tuh36651/work/Spine-Dendrite-Calcium/Mushroom_wER_wZones.ugx -numRefs 0 -ryr 1 -ryrD 3.6 -ip3r 1 -ip3rD %d -tstep 6.25e-06 -endTime 0.05 -outName /home/tuh36651/work/VietMorphSim1/Runs/Run%i -solver ILU \n',density,j);
    
    fprintf(fileID,'mpirun -np 28 ugshell -ex /home/tuh36651/work/Spine6RyR/reconstructed_spine_wER.lua -grid /home/tuh36651/work/Spine6RyR/Spine6_wER.ugx -numRefs 0 -caInflux 0.007045 -tstep 5e-6 -endTime 0.030 -outName /home/tuh36651/work/Spine6RyR/Runs/Run%i -solver GS -setting ryr -ryrDensity %d -minDef 1e-12 -freq 1',j,density);
    fclose(fileID);

end

end

