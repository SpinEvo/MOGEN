function [ang_locs_tagmode,final_encoding_matrix] = hw_fun_IterEncoding(VesLocs3D,AngTowardsCor,AngTowardsSag,off_res,number_of_iterations,modmat_phase,modmat_lookup,opts)

rng(100);

if_use_simLabEff_modulation = opts.if_use_simLabEff_modulation;
if_fft_based = opts.if_fft_based;
if_display = opts.if_display;
if_3T = opts.if_3T;

[LRUnitVec, APUnitVec] = CalcObliqueUnitVecs(AngTowardsCor,AngTowardsSag);
artery_xy_locations = transpose(Convert3DTo2DVesLocs(VesLocs3D', LRUnitVec, APUnitVec)) * 10^(-3); % in m

motion_tolerance = 4 * 10^(-3); % in m
number_of_arteries = size(artery_xy_locations,1);

% usable Hadamard order under 200
usable_Hadamard_order = [2,4,8,12,16,20,24,32,40,48,64,80,96,128,160,192];
Hadamard_order_ind = find(usable_Hadamard_order >= number_of_arteries + 1, 1);
hadamard_matrix = hadamard(usable_Hadamard_order(Hadamard_order_ind));

all_cost = zeros(1, number_of_iterations);
all_candidate_labeling_patterns = zeros(size(hadamard_matrix, 1) + 1, number_of_arteries, number_of_iterations);
labeling_pattern_cache = [];
half_spatial_cycle_cache = [];
lab_eff_cache = [];

for i = 1:number_of_iterations
    disp(['...processing potential pattern ', num2str(i), ' out of ', num2str(number_of_iterations)])
    random_order = randperm(size(hadamard_matrix, 2) - 1) + 1; 
    cur_candidate_lab_pats = hadamard_matrix(2:end, random_order(1:number_of_arteries)); % The first column does not participate in the randomization,
                                                                   % The first row is all Control and placed at the end
    cur_candidate_lab_pats = [cur_candidate_lab_pats; 
                               -1 * ones(1,number_of_arteries ); % Label
                               1 * ones(1, number_of_arteries)]; % Control
    all_candidate_labeling_patterns(:, :, i) = cur_candidate_lab_pats;

    K = false(size(cur_candidate_lab_pats, 1), 1); % logical, if be processed before
    KK = []; % indices of the matched pattern
    if size(labeling_pattern_cache, 1) > 0 % See if we have processed the labeling pattern before
        for j = 1:size(cur_candidate_lab_pats, 1)
           found_ind  = find(all(repmat(cur_candidate_lab_pats(j, :), [size(labeling_pattern_cache, 1), 1]) == labeling_pattern_cache, 2));
           if ~isempty(found_ind)
               K(j) = true;
               KK = [KK; found_ind];
           end
        end

        % for the unprocssed labeling patterns
        new_candidate_lab_pats = cur_candidate_lab_pats(~K, :);
    else
        new_candidate_lab_pats = cur_candidate_lab_pats;
    end

    half_spatial_cycle = zeros(size(cur_candidate_lab_pats, 1), 1);
    lab_eff = zeros(size(cur_candidate_lab_pats, 1), number_of_arteries);

    if any(K)
        half_spatial_cycle(K, 1) = half_spatial_cycle_cache(KK, 1);
        lab_eff(K, :) = lab_eff_cache(KK, :);
    end
    if any(~K)
        if if_fft_based % FFT-based method (OES) to find the optimal labeling parameters
            [~,new_half_spatial_cycle,new_lab_eff] = hw_fun_OES(new_candidate_lab_pats, off_res, if_use_simLabEff_modulation,...
                artery_xy_locations, modmat_phase, modmat_lookup, if_display);   
        else % NON-FFT based method (MOGEN)to find the optimal labeling parameters
            [~,new_half_spatial_cycle,new_lab_eff] = hw_fun_MOGEN(new_candidate_lab_pats, off_res, if_use_simLabEff_modulation,...
                artery_xy_locations, modmat_phase, modmat_lookup, if_display);
        end
        half_spatial_cycle(~K, 1) = new_half_spatial_cycle;                    
        lab_eff(~K, :) = new_lab_eff;

        labeling_pattern_cache = [labeling_pattern_cache; new_candidate_lab_pats];
        half_spatial_cycle_cache = [half_spatial_cycle_cache; new_half_spatial_cycle];
        lab_eff_cache = [lab_eff_cache; new_lab_eff];  
    end                  

    % now calculate the cost                        
    static_tissue = ones(size(cur_candidate_lab_pats, 1), 1);
    simulated_encoding_matrix = [lab_eff, static_tissue];
    condition = cond(simulated_encoding_matrix); 
    min_half_spatial_cycle = min(half_spatial_cycle(:));
    cost_1 = 0.5 * (1 - 1 / condition) ^ 2; 
    cost_2 = 0.5 * (1 / (min_half_spatial_cycle / motion_tolerance - 1)) ^ 2;
    all_cost(i) = cost_1 + cost_2;
end

[~, min_cost_ind] = min(all_cost);
final_encoding_matrix = all_candidate_labeling_patterns(:, :, min_cost_ind);
% Re-order: make non-selective pairs to be the first
final_encoding_matrix(1:end,:) = final_encoding_matrix([end-1:end,1:end-2],:);
% Remove the duplicate coding design
final_encoding_matrix = unique(final_encoding_matrix,'rows','stable');
clc;disp(final_encoding_matrix);

% Generate VEPCASL_OES_Setup.txt file
if if_fft_based % FFT-based method (OES) to find the optimal labeling parameters
    [ang_locs_tagmode,~,~] = hw_fun_OES(final_encoding_matrix, off_res, if_use_simLabEff_modulation,...
        artery_xy_locations, modmat_phase, modmat_lookup, false);   
else % NON-FFT based method (MOGEN)to find the optimal labeling parameters
    [ang_locs_tagmode,~,~] = hw_fun_MOGEN(final_encoding_matrix, off_res, if_use_simLabEff_modulation,...
                artery_xy_locations, modmat_phase, modmat_lookup, false);
end

if if_3T
    ang_locs_tagmode(1,:) = [0 -100 100 0]; % tag all
    ang_locs_tagmode(2,:) = [0 -100 100 1]; % control all
end

end