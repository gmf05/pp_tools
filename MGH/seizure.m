classdef seizure
% 
% seizure: Seizure object
% Part of the Point Process Toolbox by Grant Fiddyment
% 
% Constructor: sz = seizure(patient_name, seizure_name)
% The function uses the names to automatically find
% ECoG and LFP .mat files using an assumed directory structure:
%   '<patient_name>/<patient_name>_<seizure_name>_<data_type>.mat'
%   
  
  properties
    Name
    Patient
    ECoG
    LFP
    Sync
  end
  
  methods
    
    % Constructor
    function obj = seizure(patient_name, seizure_name)
      obj.Name = [patient_name ' ' seizure_name]; 
      obj.Patient = patient_name;
      
      % automatically find raw + processed files
      global DATA_DIR
      global SPIKE_THRESH_ECOG
      global SPIKE_THRESH_LFP

      OLD_DIR = pwd();
      cd([DATA_DIR '/' patient_name]);
      
      fprintf(['\nLoading seizure data (' patient_name ' ' seizure_name ')\n']);
      
      % (1) find ECoG:
      fprintf(['ECoG spike threshold: ' num2str(SPIKE_THRESH_ECOG) '\n']);
      raw_ecog_file = dir([patient_name '_' seizure_name '_ECoG.mat']);
      pp_ecog_file = dir([patient_name '_' seizure_name '_ECoG_pp_thresh' num2str(SPIKE_THRESH_ECOG) '.mat']);
      raw_ms_file = dir([patient_name '_' seizure_name '_LFP_ECoG_EEG.mat']);
      pp_lfp_file = dir([patient_name '_' seizure_name '_LFP_pp_thresh' num2str(SPIKE_THRESH_LFP) '.mat']);
      if isempty(raw_ecog_file)
        if isempty(raw_ms_file)
          if isempty(pp_ecog_file)
            error('ERROR: No ECoG data.');
          else
            fprintf('No raw ECoG data found.\n');
            ecog.RawMatFile = [];
          end
        else
          ecog.RawMatFile = raw_ms_file(1).name;
        end
      else
        ecog.RawMatFile = raw_ecog_file(1).name;
      end
      
      % (1b) if pp data is not already available, create it:
      if isempty(pp_ecog_file)
        fprintf('Point process ECoG data not found. Creating it...\n');
%         sz = mod_sz(patient_name, seizure_name);
        load(ecog.RawMatFile);        
        % create point process data object:
        data = make_seizure_pp(sz.ECoG,patient_name,seizure_name,'ECoG');
      else
        fprintf('Point process ECoG data found.\n');
        	load(pp_ecog_file(1).name);
      end
      ecog.PPData = data; 
      
      % (2) find whether multiscale data exists:      
      if isempty(raw_ms_file)
        fprintf('No raw LFP data.\n');
        lfp.RawMatFile = [];
        sync.RawMatFile = [];
      else
        lfp.RawMatFile = raw_ms_file(1).name;
        sync.RawMatFile = raw_ms_file(1).name;
      end
      
      % (2b) find whether point process LFP data exists:      
      if isempty(pp_lfp_file)
        if ~isempty(raw_ms_file)
          % make point process LFP if possible:
          fprintf('Point process LFP data not found. Creating it...\n');        
          load(lfp.RawMatFile);
          % create point process data object:
          data = make_seizure_pp(sz.LFP,patient_name,seizure_name,'LFP');
          N_channels = size(data.dn,2);
          spikes = cell(1,N_channels);
          for i = 1:N_channels
            spikes{i} = data.t(data.dn(i,:));
          end
          save([patient_name '_' seizure_name '_LFP_spikes_thresh' num2str(thresh)],'spikes');
          
        else
          data = [];
        end
      else
        fprintf('Point process LFP data found.\n');
        load(pp_lfp_file(1).name);
      end
      lfp.PPData = data;
      
      % save data
      obj.ECoG = ecog;
      obj.LFP = lfp;      
      obj.Sync = sync;
      
      % back to original directory
      cd(OLD_DIR);
      
      fprintf('Done!\n\n');
      
    end
    
    function plot(obj)
      
    
  end
  end
end
