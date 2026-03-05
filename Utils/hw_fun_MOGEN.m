function [ang_locs_tagmode,half_spatial_cycle, label_effciency] = hw_fun_MOGEN(labeling_patterns,off_res,use_simLabEff_modulation,artery_xy_locations,modmat_phase,modmat_lookup,if_display)

labeling_plane_center = [0, 0];

labeling_patterns = (labeling_patterns+1)*(pi/2);
labeling_patterns_offres = zeros(size(labeling_patterns));
% Remove the off resonance phase that will accrue between pulses due to the field inhomogeneity
for n = 1:size(labeling_patterns,1)
    labeling_patterns_offres(n,:) = labeling_patterns(n,:) - off_res;  % no offset
end
% Turn the encodings in to complex numbers (e^itheta)
labeling_patterns_offres = exp(1i*labeling_patterns_offres);

motion_tolerance = 4 * 10^(-3); % in m
number_of_arteries = size(artery_xy_locations,1);
number_of_lab_pat = size(labeling_patterns,1);

%% MOGEN Algorithm
img_size = 1024; % FOV = img_size * resolution
% CoordSysB
resolution = 1 * 10^(-3); % in m, for CoordSysB
CoordSysB_origin_in_CoordSysA = ([img_size, img_size] / 2 + 1); % x, y of CoordSysB origin in CoordSysA
artery_xy_loc = artery_xy_locations - repmat(labeling_plane_center, [number_of_arteries, 1]); 
artery_x_on_img = round(artery_xy_loc(:, 1) / resolution + CoordSysB_origin_in_CoordSysA(1));
artery_y_on_img = round(artery_xy_loc(:, 2) / resolution  + CoordSysB_origin_in_CoordSysA(2));

% define the parameter range for searching
angles_to_search = 0:1:359;
nr_angles = length(angles_to_search);
inter_artery_dist = squareform(pdist(artery_xy_loc));
diagonal_idx = eye(size(inter_artery_dist)) == 1;
min_artery_dist = min(inter_artery_dist(~diagonal_idx));
% I think that's the main reason why we generated lower spatial frequency
% compared with Tom's OES approach
step_artery_dist = 2 * 10^(-3); max_artery_dist = 150 * 10^(-3);
spatial_cycles_to_search = min_artery_dist : step_artery_dist : max_artery_dist; % in m
nr_spatial_cycles = length(spatial_cycles_to_search);

offset_range = [-75, 74] * 10^(-3); % in m
offset_step = 1 * 10^(-3);
offsets_to_search0 = (offset_range(1) - motion_tolerance) : offset_step : (offset_range(2) + motion_tolerance);
nr_offsets0 = length(offsets_to_search0);  

% calculate the projected distance (to isocenter) of the artery locations
orientation_matrix = [cos(angles_to_search * pi / 180); sin(angles_to_search * pi / 180)];
projected_dist = artery_xy_loc * orientation_matrix;

% calculate the norminal phase 
shifted_artery_loc = repmat(projected_dist, [1, 1, nr_offsets0]) + ...
    repmat(reshape(offsets_to_search0, [1, 1, nr_offsets0]), [number_of_arteries, nr_angles, 1]);
norminal_phase = repmat(shifted_artery_loc, [1, 1, 1, nr_spatial_cycles]) ./ ...
                repmat(reshape(spatial_cycles_to_search, [1, 1, 1, nr_spatial_cycles]), [number_of_arteries, nr_angles, nr_offsets0, 1])* 2 * pi - pi;

% calculate the labeling efficiency
lab_eff_to_search0 = get_labeling_efficiency(use_simLabEff_modulation,norminal_phase,modmat_phase,modmat_lookup); % Hongwei: modified for off-resonance correction

nr_steps_for_motion = round(motion_tolerance / offset_step);
H = exp(-(0:nr_steps_for_motion-1) / (nr_steps_for_motion/4));
H = [H(end:-1:2), H];
H = H / sum(H);
lab_eff_to_search0 = filter(H, 1, lab_eff_to_search0, [], 3);

K = (offsets_to_search0 >= offset_range(1)) & (offsets_to_search0 <= offset_range(2));
offsets_to_search = offsets_to_search0(K);
lab_eff_to_search = lab_eff_to_search0(:, :, nr_steps_for_motion * 2: end-1, :); % ZChen: Why not start from nr_steps_for_motion * 2 + 1?
nr_offsets = length(offsets_to_search);

label_effciency = zeros(number_of_lab_pat, number_of_arteries);
ref_xy_location = zeros(number_of_lab_pat, 2);
Gxy_angle = zeros(number_of_lab_pat, 1);
half_spatial_cycle = zeros(number_of_lab_pat, 1);
label_length = zeros(number_of_lab_pat, 1);
label_offset = zeros(number_of_lab_pat, 1);

tagang = zeros(1, number_of_lab_pat);
encs = zeros(2, number_of_lab_pat);
for i = 1 : number_of_lab_pat
    labeling_cycle = labeling_patterns_offres(i, :)';
    total_signal_to_search = sum(repmat( -1 * abs(labeling_cycle) .* (sign(real(labeling_cycle))), [1, nr_angles, nr_offsets, nr_spatial_cycles]) .* lab_eff_to_search, 1);
    max_sig_ind0 = find(total_signal_to_search(:) == max(total_signal_to_search(:)));
    [~, angle_ind0, offset_ind0, spatial_cycle_ind0] = ind2sub([1, nr_angles, nr_offsets, nr_spatial_cycles], max_sig_ind0);
    [~, tmp_ind] = max(spatial_cycles_to_search(spatial_cycle_ind0));
    angle_ind = angle_ind0(tmp_ind);
    offset_ind = offset_ind0(tmp_ind); 
    
    spatial_cycle_ind = spatial_cycle_ind0(tmp_ind);
    Gxy_angle(i) = angles_to_search(angle_ind) / 180 * pi;
    half_spatial_cycle(i) = spatial_cycles_to_search(spatial_cycle_ind) / 2;
    label_length(i) = spatial_cycles_to_search(spatial_cycle_ind) / 2;
    label_offset(i) = offsets_to_search(offset_ind) - half_spatial_cycle(i);
    label_effciency(i, :) = lab_eff_to_search(:, angle_ind, offset_ind, spatial_cycle_ind);
    ref_xy_location(i, :) = label_offset(i) * [cos(Gxy_angle(i)), sin(Gxy_angle(i))];
    encs(:,i) = [offsets_to_search(offset_ind); label_offset(i)] * 1e3; % in mm
    tagang(i) = mod(- angles_to_search(angle_ind) - 90, 360);
end

tagmode = 2 * ones(1, number_of_lab_pat);
ang_locs_tagmode = [tagang;encs;tagmode];
ang_locs_tagmode = ang_locs_tagmode';

%% Display
if if_display
    for i = 1 : number_of_lab_pat
        [xx, yy] = meshgrid(1:img_size); 
        xx = (xx - CoordSysB_origin_in_CoordSysA(1)) * resolution;
        yy = (yy - CoordSysB_origin_in_CoordSysA(2)) * resolution;

        projected_dist = xx * cos(Gxy_angle(i)) + yy *  sin(Gxy_angle(i));
        lab_eff = get_labeling_efficiency(use_simLabEff_modulation, projected_dist / half_spatial_cycle(i) * pi + label_offset(i) / half_spatial_cycle(i) * pi,modmat_phase,modmat_lookup);

        % display region
        min_ind = round(min(artery_xy_locations(:)) / resolution) * 2 + (img_size / 2 + 1); 
        max_ind = round(max(artery_xy_locations(:)) / resolution) * 2 + (img_size / 2 + 1); 

        figure, imshow(lab_eff, [-1, 1]);
        colormap(gca, jet); colorbar;
        hold on
        for j = 1 : number_of_arteries
            if labeling_patterns(i, j) == 0
                plot(artery_x_on_img(j), artery_y_on_img(j), '+k', 'markersize', 5)
            elseif labeling_patterns(i, j) == pi
                plot(artery_x_on_img(j), artery_y_on_img(j), '+w', 'markersize', 5)
            end
        end
        title(['Pattern: ', num2str(i)])
        hold off
        axis on
        xlim([min_ind, max_ind])
        ylim([min_ind, max_ind])
        set(gca, 'ydir', 'normal')
        % set(gca, 'xdir', 'reverse')
    end
end

end