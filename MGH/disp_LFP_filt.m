patient = 'MG49'; seizure = 'Seizure45';
for response = 1:10:100
  response
  open([DATA_DIR '/' patient '/' patient '_filt/' patient '_' seizure '_LFP_pp_filt_c' num2str(response) '.fig']);
  subplot(312); caxis([0,2]); ylim([0,500]); colorbar;
  subplot(313); caxis([0.9,1.1]); ylim([0,150]); colorbar;
  pause; close;

  % if LEFT, response -1
  % if RIGHT, response +1

end
