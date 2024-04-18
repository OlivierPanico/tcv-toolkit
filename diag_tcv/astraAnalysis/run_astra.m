%shot = 80257;
47778
params.kin_par.AMJ = 2;
params.kin_par.ABEAM = 2;
params.kin_par.ZBEAM = 1;
[dataD] = astra_tcv_automatic(shot, time_window, "params", params);


time2plot = [];
figure2plot = 7000;
what2plot = {'power_time', 'power_profile',  'NB_time'};
ASTRA_TCV_summary(dataH,time2plot,figure2plot,what2plot);
figure2plot = 8000;
ASTRA_TCV_summary(dataD,time2plot,figure2plot,what2plot);


itime=50;
T = dataH.out.T;
rho = dataH.out.RHOPSI;
HQIEFF = dataH.out.QIEFF;
HQEEFF = dataH.out.QEEFF;
DQIEFF = dataD.out.QIEFF;
DQEEFF = dataD.out.QEEFF;
a=HQIEFF(:,:)./(HQIEFF(:,:)+HQEEFF(:,:));
b=DQIEFF(:,:)./(DQIEFF(:,:)+DQEEFF(:,:));
figure;
plot(rho(:,itime),a(:,itime), 'b');hold on;
plot(rho(:,itime),b(:,itime), 'r');
legend('HQIEFF/(HQIEFF+HQEEFF)', 'DQIEFF/(DQIEFF+DQEEFF)');
title(sprintf('T=%d',T(itime))); 
grid;

figure;
plot(rho(:,itime),abs((a(:,itime)-b(:,itime))))


figure;
pcolor(transpose(a)-transpose(b));
