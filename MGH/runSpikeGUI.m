% load /projectnb/ecog/Data/MG49/MG49_Seizure36_LFP_ECoG_EEG.mat; tmin=50; tmax=150;
% load /projectnb/ecog/Data/MG49/MG49_Seizure45_LFP_ECoG_EEG.mat; tmin=50; tmax=150;
% load /projectnb/ecog/Data/MG63/MG63_Seizure3_LFP_ECoG_EEG.mat; tmin=0; tmax=105; 
% load /projectnb/ecog/Data/MG63/MG63_Seizure4_LFP_ECoG_EEG.mat; tmin=0; tmax=105;
load /projectnb/ecog/Data/BW09/BW09_Seizure1_LFP_ECoG_EEG.mat; tmin=0; tmax=105;
% load /projectnb/ecog/Data/BW09/BW09_Seizure2_LFP_ECoG_EEG.mat; tmin=0; tmax=105;
% load /projectnb/ecog/Data/BW09/BW09_Seizure3_LFP_ECoG_EEG.mat; % tmin=0; tmax=90;
% load /projectnb/ecog/Data/NY442/NY442_Seizure1_LFP_ECoG_EEG.mat; % tmin=0; tmax=90;
% load /projectnb/ecog/Data/RIE01/RIE01_Seizure1_LFP_ECoG_EEG.mat; % tmin=0; tmax=90;

dsz = sz.LFP.Data';
tsz = sz.LFP.Time;
t1 = getclosest(tsz,tmin);
t2 = getclosest(tsz,tmax);
T=tsz(t1:t2); D=dsz(:,t1:t2);

%%
% z-scoring:
Z=D;
for i = 1:size(D,1)
  if mod(i,30)==0, i, end
  Z(i,:)=-(Z(i,:)-mean(dsz(i,:)))./std(dsz(i,:)); % NOTE we z-score the *whole* time series dsz
end

close all;
spikeGUI(T,D,Z);