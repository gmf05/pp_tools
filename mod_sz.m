% assume sz is loaded

function sz = mod_sz(patient_name, seizure_name)
  
  global DATA_DIR
  global DO_SAVE
  
  file_name = [DATA_DIR '/' patient_name '/' patient_name '_' ...
              seizure_name '_ECoG.mat'];
  if ~exist(file_name, 'file'),
    file_name = [DATA_DIR '/' patient_name '/' patient_name '_' ...
              seizure_name '_LFP_ECoG_EEG.mat'];
  end
  load(file_name);

  if isfield(sz.ECoG, 'BadLabels')
    return;
  end
  
  bad_chan_file = [DATA_DIR '/' patient_name '/' patient_name '_old.m']; 
  if exist(bad_chan_file, 'file')
    run(bad_chan_file);
  else
    error(['ERROR: no bad labels defined in ' bad_chan_file '\n']);
  end
  
  for k = 1:length(dCell)
    if isequal(dCell(k).Seizure, seizure_name), seizure_ind = k; break; end;
  end
  
  if isfield(sz,'ECoG')
    sz.ECoG.Labels = str2cell(sz.ECoG.Labels);
    bad = dCell(seizure_ind).ECoG.BadLabels;
    badind = getnameidx(sz.ECoG.Labels, bad);
    goodind = setdiff(1:length(sz.ECoG.Labels), badind);
    sz.ECoG.Labels = {sz.ECoG.Labels{goodind}};
    sz.ECoG.BadLabels = bad;
    sz.ECoG.Data = sz.ECoG.Data(:,goodind);
    
    fprintf(['Removed ' num2str(length(badind)) ' ECoG channels\n']);
    
    sz.ECoG.Onset = sz.Onset;
    sz.ECoG.Offset = sz.Offset;
  end
  if isfield(sz,'LFP')
    sz.LFP.Labels = str2cell(sz.LFP.Labels);
    bad = dCell(seizure_ind).LFP.BadLabels;
    badind = getnameidx(sz.LFP.Labels, bad);
    goodind = setdiff(1:length(sz.LFP.Labels), badind);
    sz.LFP.Labels = {sz.LFP.Labels{goodind}};
    sz.LFP.Data = sz.LFP.Data(:,good_ind);
    
    fprintf(['Removed ' num2str(length(badind)) ' LFP channels\n']);
    
    sz.LFP.BadLabels = bad;
    sz.LFP.Onset = sz.Onset;
    sz.LFP.Offset = sz.Offset;
  end
  
  if DO_SAVE, save(file_name, 'sz'); end
  
end
