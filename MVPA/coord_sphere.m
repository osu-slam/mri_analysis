%% coord_sphere
% Creates map of searchlight spheres for each subject
% Updated to run under SLAM_isss_multi_MVPA_batch
% 
% MM/DD/YY -- CHANGELOG
% 02/14/20 -- Cloned for YA_FST, made universal

function coord_sphere(subj, study, dd)
%% Check input
if ~isstruct(subj) || length(subj) ~= 1
    error('Input (subj, study, dd) where subj is a SINGLE struct')
end

if ~isstruct(study)
    error('Input (subj, study, dd) where study is struct')
end

if ~isnumeric(dd)
    error('Input (subj, study, dd) where dd specifies which design')
end

%% Pathing and parameters
radius = study.mvpa.radius;

disp(['Radius of ' num2str(radius) ' voxels'])
dir_subj = fullfile(study.path, 'data', subj.name); 
dir_data_MVPA = fullfile(dir_subj, 'MVPA'); 

design = study.design(dd); 
dir_design = fullfile(dir_subj, 'design', design.name); 
mask_file = fullfile(dir_design, 'mask.nii'); 

%% Load mask
Vmask = spm_vol(mask_file);
mask_mat = spm_read_vols(Vmask);
mask_idx = find(mask_mat);  %%% The indices of where mask=1
num_mask_voxels = length(mask_idx);

x_size = Vmask.dim(1);
y_size = Vmask.dim(2);
z_size = Vmask.dim(3);

x_coord_vec = 1:x_size;
y_coord_vec = 1:y_size;
% z_coord_vec = 1:z_size;

%% Now make a 3D mesh-grid of x,y and z coords
% x-coord is the row in this grid
x_coord_grid = x_coord_vec' * ones(1,y_size);
% y-coord is the col in this grid
y_coord_grid = ones(x_size,1) * y_coord_vec;
% Now stack these in the 3rd dimension, to make coords cubes
x_coord_cube = [];
y_coord_cube = [];
z_coord_cube = [];

for z_slice = 1:z_size
    x_coord_cube = cat(3,x_coord_cube,x_coord_grid);
    y_coord_cube = cat(3,y_coord_cube,y_coord_grid);
    z_coord_cube = cat(3,z_coord_cube,z_slice*ones(size(x_coord_grid)));
end

% Now go through the volume and calculate sphere coords
sphere_XYZ_indices_cell = cell(num_mask_voxels,1);

% The ordering of x,y and z in zero_meaned_time_courses is this:
[x_in_mask,y_in_mask,z_in_mask]=ind2sub(size(mask_mat),mask_idx);
% XYZ has three rows, and one col for every voxel in the mask
XYZ = [x_in_mask,y_in_mask,z_in_mask]';

%% Make a predefined lookup-table of XYZ x-locations,
% to save time in the loop below
X_max = max(XYZ(1,:));
X_index_lookup_table = cell(X_max,1);
for x = 1:X_max
    X_index_lookup_table{x} = find( XYZ(1,:)==x );
end

%% Heavy lifting
parfor voxel_num =1:length(mask_idx)
    if rem(voxel_num,2000)==0
        disp([ 'Voxel number ' num2str(voxel_num) ' out of ' num2str(num_mask_voxels) ]);
    end

    [center_x, center_y, center_z] = ind2sub(size(mask_mat),mask_idx(voxel_num));

    distance_cube = sqrt( (x_coord_cube - center_x).^2 + ...
        (y_coord_cube - center_y).^2 + ...
        (z_coord_cube - center_z).^2  );

    within_sphere = ( distance_cube <= radius);
    within_sphere_and_mask = within_sphere.*mask_mat;
    sphere_inds = find(within_sphere_and_mask);
%     within_sphere_inds = find( distance_cube <= sphere_radius );
%     sphere_inds = intersect(within_sphere_inds,mask_inds);

    [sphere_x_coords, sphere_y_coords, sphere_z_coords] = ind2sub( size(mask_mat), sphere_inds);

    num_sphere_voxels = length(sphere_inds);
    XYZ_index_list = zeros(1,num_sphere_voxels);

    for sphere_vox_num = 1:num_sphere_voxels

        x = sphere_x_coords(sphere_vox_num);
        y = sphere_y_coords(sphere_vox_num);
        z = sphere_z_coords(sphere_vox_num);

        X_matches = X_index_lookup_table{x};
        XY_matches = X_matches(XYZ(2,X_matches)==y);
        XYZ_index = XY_matches(XYZ(3,XY_matches)==z);

        XYZ_index_list(sphere_vox_num) = XYZ_index;
    end   %%% End of loop through voxels within this sphere

    sphere_XYZ_indices_cell{voxel_num} = XYZ_index_list; %#ok<PFOUS>
end

%% Save file
filename = fullfile(dir_data_MVPA, [subj.name '_' design.name '_spheres_radius' num2str(study.mvpa.radius) '.mat']); 
save(filename, 'sphere_XYZ_indices_cell');

end
