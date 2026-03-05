function [ang_locs_tagmode,half_spatial_cycle, label_effciency] = hw_fun_OES(labeling_patterns,off_res,use_simLabEff_modulation,artery_xy_locations,modmat_phase,modmat_lookup,if_display)

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
minimum_spatial_cycle = motion_tolerance * 4; % in m
number_of_arteries = size(artery_xy_locations,1);
number_of_lab_pat = size(labeling_patterns,1);

%% OES Algorithm
% CoordSysA
img_size = 1024; % FOV = img_size * resolution
img = zeros(img_size, img_size);
% CoordSysB
resolution = 1 * 10^(-3); % in m, for CoordSysB
CoordSysB_origin_in_CoordSysA = ([img_size, img_size] / 2 + 1); % x, y of CoordSysB origin in CoordSysA
freq_resolution = 1 / (img_size * resolution); % in m^(-1), frequency (i.e. kspace) resolution
KCoordSysB_origin_in_KCoordSysA = ([img_size, img_size] / 2 + 1); % kx, ky of KCoordSysB in KCoordSysA

min_spatial_cycle = minimum_spatial_cycle; % in m
max_radius_I = ceil((1 / min_spatial_cycle) / freq_resolution);
[JJ, II] = meshgrid(1:img_size);
position_radius = sqrt((II - KCoordSysB_origin_in_KCoordSysA(2)).^2 + (JJ - KCoordSysB_origin_in_KCoordSysA(1)).^2);
position_weighting = (1 - position_radius / max(max(position_radius))) / 2;

% fft image mask
% Hongwei: The centre point has been excluded, where position_radius == 0
img_fft_mask = position_radius ~= 0 & position_radius < 100; 
img_fft_mask(position_radius > max_radius_I) = 0;

position_weighting(~img_fft_mask) = 0;

offset_artery_xy_locations = artery_xy_locations - repmat(labeling_plane_center, [number_of_arteries, 1]); 
artery_x_on_img = round(offset_artery_xy_locations(:, 1) / resolution + CoordSysB_origin_in_CoordSysA(1));
artery_y_on_img = round(offset_artery_xy_locations(:, 2) / resolution  + CoordSysB_origin_in_CoordSysA(2));

artery_indice_on_img = sub2ind(size(img), artery_y_on_img, artery_x_on_img);
encs = zeros(2, number_of_lab_pat);
tagang = zeros(1, number_of_lab_pat);
max_freq_phase = zeros(1, number_of_lab_pat);
half_spatial_cycle = zeros(number_of_lab_pat,1);
offset_dist = zeros(1, number_of_lab_pat);
max_kspace_val_kxy = zeros(number_of_lab_pat,2);
label_effciency = zeros(number_of_lab_pat, number_of_arteries);

for i = 1 : number_of_lab_pat
    for j = 1 : number_of_arteries
        img(artery_indice_on_img(j)) = labeling_patterns_offres(i, j);
    end
    
    img_fft = fftshift(fft2(fftshift(img)));
    
    img_fft(~img_fft_mask) = 0; % apply mask
    img_fft_mag = abs(img_fft);
    img_fft_mag = img_fft_mag / max(img_fft_mag(:));
    [~, max_mag_ind2] = max(img_fft_mag(:) + position_weighting(:));
    [max_mag_I2, max_mag_J2] = ind2sub(size(img_fft), max_mag_ind2);
    max_kspace_val = img_fft(max_mag_I2, max_mag_J2); % May have some error caused by smoothing
    max_freq_phase(i) = angle(max_kspace_val) + pi;

    max_kspace_val_kxy(i,1) = (max_mag_J2 - KCoordSysB_origin_in_KCoordSysA(1)) * freq_resolution;
    max_kspace_val_kxy(i,2) = (max_mag_I2 - KCoordSysB_origin_in_KCoordSysA(2)) * freq_resolution;
    dist_to_freq_center = norm(max_kspace_val_kxy(i,:));
    spatial_cycle = 1 / dist_to_freq_center; % in m
    half_spatial_cycle(i) = spatial_cycle / 2;
    offset_dist(i) = - max_freq_phase(i) / (2 * pi) * spatial_cycle * 1e3; % in mm
    
    tagang(i) = atan2(max_kspace_val_kxy(i,1),max_kspace_val_kxy(i,2)) * 180 / pi;
    encs(:,i) = [offset_dist(i) - half_spatial_cycle(i) * 1e3; offset_dist(i)]; % in mm
    
    % calculate the labeling efficiency
    for j = 1 : number_of_arteries
        cur_x = (artery_x_on_img(j) - CoordSysB_origin_in_CoordSysA(1)) * resolution;
        cur_y = (artery_y_on_img(j) - CoordSysB_origin_in_CoordSysA(2)) * resolution;
        cur_projected_dist = (cur_x * max_kspace_val_kxy(i,1) + cur_y * max_kspace_val_kxy(i,2)) / norm(max_kspace_val_kxy(i,:));
        label_effciency(i, j) = get_labeling_efficiency(use_simLabEff_modulation,cur_projected_dist / spatial_cycle * 2 * pi + max_freq_phase(i),modmat_phase,modmat_lookup);
    end
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

        projected_dist = (xx * max_kspace_val_kxy(i,1) + yy * max_kspace_val_kxy(i,2)) / norm(max_kspace_val_kxy(i,:));

        lab_eff = get_labeling_efficiency(use_simLabEff_modulation, projected_dist / half_spatial_cycle(i) * pi + max_freq_phase(i),modmat_phase,modmat_lookup);

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