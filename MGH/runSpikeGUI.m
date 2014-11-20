% load /projectnb/ecog/Data/MG49/MG49_Seizure36_LFP_ECoG_EEG.mat
% load /projectnb/ecog/Data/MG49/MG49_Seizure45_LFP_ECoG_EEG.mat % 
load /projectnb/ecog/Data/MG63/MG63_Seizure3_LFP_ECoG_EEG.mat % ~30-103
% load /projectnb/ecog/Data/MG63/MG63_Seizure4_LFP_ECoG_EEG.mat % ~30-100?
% load /projectnb/ecog/Data/BW09/BW09_Seizure1_LFP_ECoG_EEG.mat
% load /projectnb/ecog/Data/BW09/BW09_Seizure2_LFP_ECoG_EEG.mat
% load /projectnb/ecog/Data/BW09/BW09_Seizure3_LFP_ECoG_EEG.mat
dsz = sz.LFP.Data';
tsz = sz.LFP.Time;

% %%
% figure, plot(t,d(2,:))

%%
% tmin=50; tmax=150; % mg49
% tmin=110; tmax=145; %mg49-s36
% tmin=115; tmax=150; %mg49-s45
tmin = 0; tmax = 100;
% tmin=50; tmax=103; % mg63-s3
% tmin=30; tmax=105; % mg63-s4
% tmin=20; tmax=90; % bw09-s3
% tmin=10; tmax=100; % bw09-s4
% tmin=70; tmax=140;

t1 = getclosest(tsz,tmin);
t2 = getclosest(tsz,tmax);
T=tsz(t1:t2); D=dsz(:,t1:t2);
% % z-scoring:
Z=D;
for i = 1:size(D,1), if mod(i,30)==0, i, end, Z(i,:)=-(Z(i,:)-mean(dsz(i,:)))./std(dsz(i,:)); end

close all;
spikeGUI(T,D,Z);