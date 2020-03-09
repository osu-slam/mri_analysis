%% Convert_FM_VDM
% Converts fieldmap (from FSL) to voxel displacement map (for SPM)
%
% MM/DD/YY: CHANGELOG
% 02.10.20: File initialized. 

% Make changes to IP
% Plug it back in

function convert_fm_vdm(dir_subj, study)

load FieldMap_IP.mat

IP.pP.fname = fullfile(dir_subj, 'realign', 'fpm_fieldmap.nii'); 
V = spm_vol(IP.pP.fname); 
IP.pP = V; 

IP.fm.fpm = spm_read_vols(IP.pP);
IP.fm.jac = pm_diff(IP.fm.fpm,2);
IP.tert   = study.scan.tert; 

FieldMap('fm2vdm', IP)

end