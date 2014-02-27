pt_name = 'MG49';  sz_name = 'Seizure45';
load([pt_name '_' sz_name '_pp_thresh1_static.mat']);
sz = seizure(pt_name,sz_name);
d = sz.LFP.PPData;
d = d.downsample(32);
XI_coarse = xi(ms,d);
save([pt_name '_' sz_name '_XI_coarse'],'XI_coarse','p');
