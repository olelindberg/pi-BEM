clear all;close all;
%% KVLCC2 PMM app-shallow FHR

% 			GLOBAL POSITION			VELOCITY			ACCELERATION								PROPELLER				RUDDER									
% name	type	flag	x	y	psi	u	v	r	u'	v'	r'	X	Y	N	SINKAGE	trim	given	measured	T	Q	given		measured		X	Y	Q	K	Drift angle	Uc
%			(m)	(m)	(°)	(m/s)	(m/s)	(°/s)	(m/s2)	(m/s2)	(°/s2)	(N)	(N)	(Nm)	(mm)	(mm/m)	(rpm)	(rpm)	(N)	(Nmm)	(°)		(°)		(N)	(N)	(Nmm)	(Nm)	(°)	(m/s)

%% Sheet 1 - depth = 0.499
[data, run_ID, alldata] = xlsread('App08_KVLCC_PMM_results.xls', 3, 'A146:AB747');
PMMY2 = data(1:48,2:end);       % Harmonic sway PMMY2
PMMPSI2 = data(50:602,2:end);   % Harmonic yaw PMMPSI2
save('Data_KVLCC2_FHR_PMM_Shallow_499mm.mat','PMMY2','PMMPSI2');

clear all;
%% Sheet 2 - depth = 0.416
[data, run_ID, alldata] = xlsread('App08_KVLCC_PMM_results.xls', 4, 'A146:AB747');
PMMY2 = data(1:48,2:end);       % Harmonic sway PMMY2
PMMPSI2 = data(50:602,2:end);   % Harmonic yaw PMMPSI2
save('Data_KVLCC2_FHR_PMM_Shallow_416mm.mat','PMMY2','PMMPSI2');

clear all;
%% Sheet 3 - depth = 0.333
[data, run_ID, alldata] = xlsread('App08_KVLCC_PMM_results.xls', 5, 'A146:AB747');
PMMY2 = data(1:48,2:end);       % Harmonic sway PMMY2
PMMPSI2 = data(50:602,2:end);   % Harmonic yaw PMMPSI2
save('Data_KVLCC2_FHR_PMM_Shallow_333mm.mat','PMMY2','PMMPSI2');