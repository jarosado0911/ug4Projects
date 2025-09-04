function calciumDataPlot()

fig=figure('units','normalized','outerposition',[0 0 1 1]);
%set(fig,'Visible', 'off');
v = VideoWriter('Spine11RyRSercaNeck1_sweep','MPEG-4');
open(v)
for j = 0:500
    clf
fileName1 = sprintf('Runs/Run%i/meas/data_meas_dend_ca_cyt',j);
fileName2 = sprintf('Runs/Run%i/meas/data_meas_neck_ca_cyt',j);
fileName3 = sprintf('Runs/Run%i/meas/data_meas_head_ca_cyt',j);
fileName4 = sprintf('Runs/Run%i/meas/data_er_ca_er',j);
fileID1 = fopen(fileName1,'r'); fileID2 = fopen(fileName2,'r');
fileID3 = fopen(fileName3,'r'); fileID4 = fopen(fileName4,'r');

sizeA = [2 Inf]; sizeB = [2 Inf]; sizeC = [2 Inf]; sizeD = [2 Inf];

A=fscanf(fileID1,'%f %f', sizeA); B=fscanf(fileID2,'%f %f', sizeB);
C=fscanf(fileID3,'%f %f', sizeC); D=fscanf(fileID4,'%f %f', sizeD);
fprintf('Read in file number %i\n',j);
fclose('all');

%if A(1,end)<=0.015
%    continue
%end
    density = 0 + j*0.01;
    subplot(1,2,1)
    hold on
    plot(C(1,:),C(2,:),'color',[0 0.5 0],'LineWidth',2);
    plot(B(1,:),B(2,:),'b','LineWidth',2);
    plot(A(1,:),A(2,:),'r','LineWidth',2);
    hold off
    xlim([0 0.03])
    ylim([0 1e-5])
    legend('Head','Neck','Dend')
    xlabel('Time [seconds]')
    ylabel('[Ca^{2+}] mol/l')
    title(sprintf(' caInflux = 4.119e-18 mol/s.um^2'));
    set(gca, 'FontSize', 16)
    
    subplot(1,2,2)
    plot(D(1,:),D(2,:),'b','LineWidth',2);
    xlim([0 0.03])
    ylim([0 3e-4])
    xlabel('Time [seconds]')
    ylabel('[Ca^{2+}] mol/l')
    title('ER Calcium Conc.')
    set(gca, 'FontSize', 16)
    
    sgtitle(sprintf('Calcium Concentration RyR density = %0.3f um^{-2}',density),'fontsize',18);
    drawnow
    thisframe=getframe(fig);
    writeVideo(v, thisframe);
end

end

