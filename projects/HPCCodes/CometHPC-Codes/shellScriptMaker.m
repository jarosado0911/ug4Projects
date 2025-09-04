function shellScriptMaker(nFiles)
nn=linspace(0,nFiles,nFiles+1);
density = -0.01;

for k=1:length(nn)
    j=k-1;
    filename=sprintf('Runs/%s_%i.sb','script',j);
    fileID = fopen(filename,'wt');
    density = density + 0.01;
    fprintf(fileID,'#!/bin/sh \n');
    fprintf(fileID,sprintf('#SBATCH --job-name="S11-Neck1-%i" \n',k));
    fprintf(fileID,'#SBATCH --output="Spine11-Neck2.out" \n');
    fprintf(fileID,'#SBATCH --partition=compute \n');
    fprintf(fileID,'#SBATCH --nodes=1 \n');
    fprintf(fileID,'#SBATCH --ntasks-per-node=24 \n');
    fprintf(fileID,'#SBATCH --export=ALL \n');
    fprintf(fileID,'#SBATCH -t 01:20:00 \n\n');
    
    %fprintf(fileID,'cd $PBS_O_WORKDIR \n\n');
    
    %fprintf(fileID,'mpirun -np 28 ugshell -ex /home/tuh36651/work/Spine-Dendrite-Calcium/reconstructed_spine_wER.lua -grid /home/tuh36651/work/Spine-Dendrite-Calcium/Mushroom_wER_wZones.ugx -numRefs 0 -ryr 1 -ryrD 3.6 -ip3r 1 -ip3rD %d -tstep 6.25e-06 -endTime 0.05 -outName /home/tuh36651/work/VietMorphSim1/Runs/Run%i -solver ILU \n',density,j);
    
    fprintf(fileID,'ibrun -np 24 ugshell -ex /home/tuh36651/ug4/runs/Spine11ryrNeck1/reconstructed_spine_wER.lua -grid /home/tuh36651/ug4/runs/Spine11ryrNeck1/Spine11_wER_v1.ugx -numRefs 0 -caInflux 0.0041195 -tstep 5e-6 -endTime 0.030 -outName /home/tuh36651/ug4/runs/Spine11ryrNeck1/Runs/Run%i -solver GS -setting ryr -ryrDensity %d -minDef 1e-12 -freq 1',j,density);
    fclose(fileID);

end

end

