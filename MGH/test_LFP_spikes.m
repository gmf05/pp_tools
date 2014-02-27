
% load data
patient_name = 'MG49';
seizure_name = 'Seizure45';
filename = [DATA_DIR '/' patient_name '/' patient_name ...
  '_' seizure_name '_LFP_ECoG_EEG.mat'];
load(filename);

%% pre-process LFP
t = sz.LFP.Time;
response = 1;
dpre = sz.LFP.Data(:,response);
dpost = preprocessing(dpre,t,'LFP');

subplot(211); plot(t,dpre);
subplot(212); plot(t,dpost);

%% pre-process ECoG
t = sz.ECoG.Time;
response = 1;
dpre = sz.ECoG.Data(:,response);
dpost = preprocessing(dpre,t,'LFP');

subplot(211); plot(t,dpre);
subplot(212); plot(t,dpost);


%% compare spectra before/after

run('~/Code/matlab/get_gmf_env'); % add chronux

% N_tapers = 1;
N_tapers = 5;
hi_fq = 250;

% set chronux parameters
% c_params.tapers = [hi_fq window_size/dt 2*hi_fq*window_size-N_tapers];
c_params.tapers = [5 9];
c_params.pad = 0;
c_params.Fs = 1/(t(2)-t(1));
c_params.fpass = [0 hi_fq];
c_params.err = 0;
c_params.trialave = [];

[Spre,fpre] = mtspectrumc(dpre,c_params);
[Spost,fpost] = mtspectrumc(dpost,c_params);
subplot(211); plot(fpre, 10*log10(Spre));
subplot(212); plot(fpost, 10*log10(Spost));

%% find spikes
clf;
plot(t,dpost); hold on;

threshs = [0.5, 1, 2];
styles = {'ro','gs','mx'};
for i=1:3
  spk_ind = hilbertspike(dpost,threshs(i),MIN_REFRACT);
  plot(t(spk_ind),dpost(spk_ind),styles{i});
  pause;
end

