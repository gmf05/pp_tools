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
      OLD_DIR = pwd();
      cd([DATA_DIR '/' patient_name]);
      
      fprintf(['\nLoading seizure data (' patient_name ' ' seizure_name ')\n']);
      
      % get ECoG:
      dECoG = get_spikes(patient_name,seizure_name,'ECoG');
      ecog.PPData = dECoG; 
      
      
      % get LFP(?):
      dLFP = get_spikes(patient_name,seizure_name,'LFP');
      lfp.PPData = dLFP;
      
      % get MUA(?):
%       dMUA = get_spikes(patient_name,seizure_name,'MUA');
%       mua.PPData = dMUA;
      
      % get EEG(?):
%       dEEG = get_spikes(patient_name,seizure_name,'EEG');
%       eeg.PPData = dEEG;
      
      % save data
      obj.ECoG = ecog;
      obj.LFP = lfp;      
%       obj.MUA = mua;
%       obj.EEG = eeg;
      obj.Sync = sync;
      
      cd(OLD_DIR);
      fprintf('Done!\n\n');
    end
    
    function plot(obj)
      
    
    end
  end
end
