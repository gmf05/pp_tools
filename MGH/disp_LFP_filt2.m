patient = 'MG49'; seizure = 'Seizure45';
for response = 1:96
  response
  open([patient '_' seizure '_LFP_pp_filt2_c' num2str(response) '.fig']);
  subplot(412); caxis([0,3]); ylim([0,500]); colorbar;
  subplot(413); caxis([0.9,1.1]); ylim([0,150]); colorbar;
  subplot(414); caxis([0.9,1.1]); ylim([0,150]); colorbar;
  pause; close;

  % if LEFT, response -1
  % if RIGHT, response +1

end
