%Estimate flux partition in mix nbi 1/2
shot = 80328
time_window=[0.6:0.01:2]
[dataD] = astra_tcv_automatic(shot, time_window);%, "params", params);
time2plot = [];
%figure2plot = 7000;
%what2plot = {'power_time', 'power_profile',  'NB_time'};
%ASTRA_TCV_summary(dataH,time2plot,figure2plot,what2plot);

shot = 80328
time_window=[0.6:0.01:2]
params=astra_tcv_params(); %changing parameter to hydrogen
params.kin_par.AMJ = 1;
params.kin_par.ABEAM = 1;
params.kin_par.ZBEAM = 1;
[dataH] = astra_tcv_automatic(shot, time_window, "params", params);
time2plot = [];
%figure2plot = 7000;
%what2plot = {'power_time', 'power_profile',  'NB_time'};
%ASTRA_TCV_summary(dataH,time2plot,figure2plot,what2plot);
TD = dataD.out.T;
rhoD = dataD.out.RHOPSI;
DQIEFF = dataD.out.QIEFF;
DQEEFF = dataD.out.QEEFF;
Dflux_partition=DQIEFF(:,:)./(DQIEFF(:,:)+DQEEFF(:,:));

TH = dataH.out.T;
rhoH = dataH.out.RHOPSI;
HQIEFF = dataH.out.QIEFF;
HQEEFF = dataH.out.QEEFF;
Hflux_partition=HQIEFF(:,:)./(HQIEFF(:,:)+HQEEFF(:,:));

itime1 = 61;
itime2 = 127;
figure;

plot(TD,Dflux_partition(43,:), 'b');hold on;
plot(TH,Hflux_partition(43,:), 'r');hold on;
%plot(,itime1),a(43,:), 'b');
legend('deuterium', 'hydrogen');
ylabel('QIEFF/(QIEFF+QEEFF)');
title(sprintf('rho=%d',mean(rhoD(43,:)))); 
grid;




%comparison H/D fluxes

shot = 81069
time_window=[0.6:0.01:1]
[dataD] = astra_tcv_automatic(shot, time_window);%, "params", params);
time2plot = [];
%figure2plot = 7000;
%what2plot = {'power_time', 'power_profile',  'NB_time'};
%ASTRA_TCV_summary(dataH,time2plot,figure2plot,what2plot);

shot = 80949
time_window=[0.6:0.01:1]
params=astra_tcv_params(); %changing parameter to hydrogen
% params.kin_par.AMJ = 1;
% params.kin_par.ABEAM = 1;
% params.kin_par.ZBEAM = 1;
[dataH] = astra_tcv_automatic(shot, time_window);%, "params", params);
time2plot = [];
%figure2plot = 7000;
%what2plot = {'power_time', 'power_profile',  'NB_time'};
%ASTRA_TCV_summary(dataH,time2plot,figure2plot,what2plot);
TD = dataD.out.T;
rhoD = dataD.out.RHOPSI;
DQIEFF = dataD.out.QIEFF;
DQEEFF = dataD.out.QEEFF;
Dflux_partition=DQIEFF(:,:)./(DQIEFF(:,:)+DQEEFF(:,:));

TH = dataH.out.T;
rhoH = dataH.out.RHOPSI;
HQIEFF = dataH.out.QIEFF;
HQEEFF = dataH.out.QEEFF;
Hflux_partition=HQIEFF(:,:)./(HQIEFF(:,:)+HQEEFF(:,:));


figure;

plot(TD,Dflux_partition(43,:), 'b');hold on;
plot(TH,Hflux_partition(43,:), 'r');hold on;
%plot(,itime1),a(43,:), 'b');
legend('deuterium', 'hydrogen');
ylabel('QIEFF/(QIEFF+QEEFF)');
title(sprintf('rho=%d',mean(rhoD(43,:)))); 
grid;


figure;
irho=50
plot(TD,Dflux_partition(irho,:), 'b');hold on;
plot(TH,Hflux_partition(irho,:), 'r');hold on;
%plot(,itime1),a(43,:), 'b');
legend('deuterium', 'hydrogen');
ylabel('QIEFF/(QIEFF+QEEFF)');
title(sprintf('rho=%d',mean(rhoD(irho,:)))); 
grid;




%comparison H/D fluxes

shot = 81084
time_window=[0.6:0.01:1]
[dataD] = astra_tcv_automatic(shot, time_window);%, "params", params);
time2plot = [];
%figure2plot = 7000;
%what2plot = {'power_time', 'power_profile',  'NB_time'};
%ASTRA_TCV_summary(dataH,time2plot,figure2plot,what2plot);

shot = 80940
time_window=[0.6:0.01:1]
params=astra_tcv_params(); %changing parameter to hydrogen
params.kin_par.AMJ = 1;
% params.kin_par.ABEAM = 1;
% params.kin_par.ZBEAM = 1;
[dataH] = astra_tcv_automatic(shot, time_window);%, "params", params);
time2plot = [];
%figure2plot = 7000;
%what2plot = {'power_time', 'power_profile',  'NB_time'};
%ASTRA_TCV_summary(dataH,time2plot,figure2plot,what2plot);
TD = dataD.out.T;
rhoD = dataD.out.RHOPSI;
DQIEFF = dataD.out.QIEFF;
DQEEFF = dataD.out.QEEFF;
Dflux_partition=DQIEFF(:,:)./(DQIEFF(:,:)+DQEEFF(:,:));

TH = dataH.out.T;
rhoH = dataH.out.RHOPSI;
HQIEFF = dataH.out.QIEFF;
HQEEFF = dataH.out.QEEFF;
Hflux_partition=HQIEFF(:,:)./(HQIEFF(:,:)+HQEEFF(:,:));


figure;
irho=43
plot(TD,Dflux_partition(irho,:), 'b');hold on;
plot(TH,Hflux_partition(irho,:), 'r');hold on;
%plot(,itime1),a(43,:), 'b');
legend('deuterium', 'hydrogen');
ylabel('QIEFF/(QIEFF+QEEFF)');
title(sprintf('rho=%d',mean(rhoD(irho,:)))); 
grid;


figure;
irho=10
plot(TD,Dflux_partition(irho,:), 'b');hold on;
plot(TH,Hflux_partition(irho,:), 'r');hold on;
%plot(,itime1),a(43,:), 'b');
legend('deuterium', 'hydrogen');
ylabel('QIEFF/(QIEFF+QEEFF)');
title(sprintf('rho=%d',mean(rhoD(irho,:)))); 
grid;


figure;
irho=40;
plot(TD,DQIEFF(irho,:), 'r');hold on;
plot(TD,DQEEFF(irho,:), 'b');hold on;
plot(TH,HQIEFF(irho,:),'r+:');hold on;
plot(TH,HQEEFF(irho,:),'b+:');hold on;
%plot(TD,Dflux_partition(irho,:), 'k');hold on;
legend('ion', 'electron', 'ion/tot');
ylabel('QIEFF/(QIEFF+QEEFF)');
title(sprintf('rho=%d',mean(rhoD(irho,:)))); 
grid;



figure;
itime=50;
plot(rhoD(:,itime),(DQIEFF(:,itime)), 'r');hold on;
plot(rhoD(:,itime),(DQEEFF(:,itime)), 'b');hold on;
plot(rhoH(:,itime),(HQIEFF(:,itime)),'r+:');hold on;
plot(rhoH(:,itime),(HQEEFF(:,itime)),'b+:');hold on;
legend('ion deuterium', 'electron deuterium', 'ion hydrogen', 'electron hydrogen');
ylabel('fluxes');
xlabel('rho');
title(sprintf('time=%d',(TD(itime)))); 
grid;



%comparison H/D fluxes

shot = 81084
time_window=[1:0.01:2]
[dataD] = astra_tcv_automatic(shot, time_window);%, "params", params);
time2plot = [];
%figure2plot = 7000;
%what2plot = {'power_time', 'power_profile',  'NB_time'};
%ASTRA_TCV_summary(dataH,time2plot,figure2plot,what2plot);

shot = 80949
time_window=[1:0.01:2]
params=astra_tcv_params(); %changing parameter to hydrogen
params.kin_par.AMJ = 1;
% params.kin_par.ABEAM = 1;
% params.kin_par.ZBEAM = 1;
[dataH] = astra_tcv_automatic(shot, time_window);%, "params", params);
time2plot = [];
%figure2plot = 7000;
%what2plot = {'power_time', 'power_profile',  'NB_time'};
%ASTRA_TCV_summary(dataH,time2plot,figure2plot,what2plot);
TD = dataD.out.T;
rhoD = dataD.out.RHOPSI;
DQIEFF = dataD.out.QIEFF;
DQEEFF = dataD.out.QEEFF;
Dflux_partition=DQIEFF(:,:)./(DQIEFF(:,:)+DQEEFF(:,:));
 
TH = dataH.out.T;
rhoH = dataH.out.RHOPSI;
HQIEFF = dataH.out.QIEFF;
HQEEFF = dataH.out.QEEFF;
Hflux_partition=HQIEFF(:,:)./(HQIEFF(:,:)+HQEEFF(:,:));


figure;
irho=43
plot(TD,Dflux_partition(irho,:), 'b');hold on;
plot(TH,Hflux_partition(irho,:), 'r');hold on;
%plot(,itime1),a(43,:), 'b');
legend('deuterium', 'hydrogen');
ylabel('QIEFF/(QIEFF+QEEFF)');
title(sprintf('rho=%d',mean(rhoD(irho,:)))); 
grid;


figure;
irho=10
plot(TD,Dflux_partition(irho,:), 'b');hold on;
plot(TH,Hflux_partition(irho,:), 'r');hold on;
%plot(,itime1),a(43,:), 'b');
legend('deuterium', 'hydrogen');
ylabel('QIEFF/(QIEFF+QEEFF)');
title(sprintf('rho=%d',mean(rhoD(irho,:)))); 
grid;


figure;
irho=40;
plot(TD,DQIEFF(irho,:), 'r');hold on;
plot(TD,DQEEFF(irho,:), 'b');hold on;
plot(TH,HQIEFF(irho,:),'r+:');hold on;
plot(TH,HQEEFF(irho,:),'b+:');hold on;
%plot(TD,Dflux_partition(irho,:), 'k');hold on;
legend('ion', 'electron', 'ion/tot');
ylabel('QIEFF/(QIEFF+QEEFF)');
title(sprintf('rho=%d',mean(rhoD(irho,:)))); 
grid;



figure;
itime=50;
plot(rhoD(:,itime),(DQIEFF(:,itime)), 'r');hold on;
plot(rhoD(:,itime),(DQEEFF(:,itime)), 'b');hold on;
plot(rhoH(:,itime),(HQIEFF(:,itime)),'r+:');hold on;
plot(rhoH(:,itime),(HQEEFF(:,itime)),'b+:');hold on;
legend('ion deuterium', 'electron deuterium', 'ion hydrogen', 'electron hydrogen');
ylabel('fluxes');
xlabel('rho');
title(sprintf('time=%d',(TD(itime)))); 
grid;

