%shot = 80257;
%47778
params=astra_tcv_params();
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
plot(,itime1),a(43,:), 'b');
legend('deuterium', 'hydrogen');
ylabel('QIEFF/(QIEFF+QEEFF)');
title(sprintf('rho=%d',mean(rho(43,:)))); 
grid;


%Comparison fluxes D/H for matching profiles
shot = 80257;
time_window=[1.55:0.01:1.65];
[dataD] = astra_tcv_automatic(shot, time_window);%, "params", params);

shot = 80931;
time_window=[0.65:0.01:0.85];
[dataH] = astra_tcv_automatic(shot, time_window);%, "params", params);

itime=50;
TH = dataH.out.T;
rhoH = dataH.out.RHOPSI;
TD = dataD.out.T;
rhoD = dataD.out.RHOPSI;
HQIEFF = dataH.out.QIEFF;
HQEEFF = dataH.out.QEEFF;
DQIEFF = dataD.out.QIEFF;
DQEEFF = dataD.out.QEEFF;
a=HQIEFF(:,:)./(HQIEFF(:,:)+HQEEFF(:,:));
b=DQIEFF(:,:)./(DQIEFF(:,:)+DQEEFF(:,:));

figure;
itimeD = 30
itimeH = 35
plot(rhoH(:,itimeH),a(:,itimeH), 'r');hold on;
plot(rhoD(:,itimeD),b(:,itimeD), 'b');
legend('HQIEFF/(HQIEFF+HQEEFF)', 'DQIEFF/(DQIEFF+DQEEFF)');
title(sprintf('T=%d, %d',TD(itimeD), TH(itimeH))); 
grid;

fig



shot = 80257;
time_window=[0.6:0.1:1.2];
[dataD] = astra_tcv_automatic(shot, time_window);%, "params", params);

shot = 80946;
time_window=[0.65:0.01:0.85];
[dataH] = astra_tcv_automatic(shot, time_window);%, "params", params);

itime=50;
TH = dataH.out.T;
rhoH = dataH.out.RHOPSI;
TD = dataD.out.T;
rhoD = dataD.out.RHOPSI;
HQIEFF = dataH.out.QIEFF;
HQEEFF = dataH.out.QEEFF;
DQIEFF = dataD.out.QIEFF;
DQEEFF = dataD.out.QEEFF;
a=HQIEFF(:,:)./(HQIEFF(:,:)+HQEEFF(:,:));
b=DQIEFF(:,:)./(DQIEFF(:,:)+DQEEFF(:,:));

figure;
itimeD = 8;
itimeH = 36;
plot(rhoH(:,itimeH),a(:,itimeH), 'r');hold on;
plot(rhoD(:,itimeD),b(:,itimeD), 'b');
legend('HQIEFF/(HQIEFF+HQEEFF)', 'DQIEFF/(DQIEFF+DQEEFF)');
title(sprintf('T=%d, %d',TD(itimeD), TH(itimeH))); 
grid;


shot = 80322;
time_window=[0.6:0.01:1.2];
[dataD] = astra_tcv_automatic(shot, time_window);%, "params", params);

shot = 80946;
time_window=[0.65:0.01:0.85];
[dataH] = astra_tcv_automatic(shot, time_window);%, "params", params);

TH = dataH.out.T;
rhoH = dataH.out.RHOPSI;
TD = dataD.out.T;
rhoD = dataD.out.RHOPSI;
HQIEFF = dataH.out.QIEFF;
HQEEFF = dataH.out.QEEFF;
DQIEFF = dataD.out.QIEFF;
DQEEFF = dataD.out.QEEFF;
a=HQIEFF(:,:)./(HQIEFF(:,:)+HQEEFF(:,:));
b=DQIEFF(:,:)./(DQIEFF(:,:)+DQEEFF(:,:));

figure;
itimeD = 61;
itimeH = 36;
plot(rhoH(:,itimeH),a(:,itimeH), 'r');hold on;
plot(rhoD(:,itimeD),b(:,itimeD), 'b');
legend('HQIEFF/(HQIEFF+HQEEFF)', 'DQIEFF/(DQIEFF+DQEEFF)');
title(sprintf('T=%d, %d',TD(itimeD), TH(itimeH))); 
grid;
