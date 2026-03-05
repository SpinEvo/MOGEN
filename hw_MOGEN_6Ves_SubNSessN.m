clc;clear;close all;
workingDir = 'E:\Hongwei\***'; cd(workingDir);
util_path = [workingDir,filesep,'Utils']; addpath(util_path);

if_fft_based = false; % true - IOES; false - MOGEN
Hadamard_order = 12; % 4,8,12,16,20,24,32,40,48,64,...
%% Arteries Siemens Coordinates
%##########################################################################
AngTowardsCor = 0; AngTowardsSag = 0; off_res = 0; isUnipolar = 0;

% L-positive, R-negative
% A-positive, P-negative
% H-positive, F-negative

VesLocs3DStr = [{'R28 , A29'}       % RECA
                {'R27 , A21'}       % RICA
                {'R15 , A14'}        % RVA
                {'L24 , A34'}       % LECA
                {'L33 , A25'}        % LICA
                {'L15 , A21'}];      % LVA

% PLEASE do change the labeling plane: [Geometry -> Saturation(2)]
z_value = hw_convertCoordinateValue('F72');

VesLocs3D = hw_convertAllCoordinates(VesLocs3DStr);

% Hongwei: reorder-RICA,RVA,LICA,LVA,RECA,LECA
VesLocs3D = VesLocs3D([2,3,5,6,1,4],:);
VesLocs3D = [VesLocs3D, z_value * ones(size(VesLocs3D, 1), 1)];

%##########################################################################


%% Preparing for 6Ves Encodings 
% Please do not change anything!!!

load([util_path,filesep,'modmat_phase.mat']);
vendor = 'Para2'; % 'Para1','Para2'

if strcmp(vendor,'Para2')
    load([util_path,filesep,'modmat_lookup_Para2_3T_bipolar.mat']);
elseif strcmp(vendor,'Para1')
    load([util_path,filesep,'modmat_lookup_Para1_3T_bipolar.mat']);
end


number_of_iterations = 1000;
opts.if_use_simLabEff_modulation = true;
opts.if_fft_based = if_fft_based;
opts.if_display = false;
opts.if_3T = true;
opts.Hadamard_order = Hadamard_order;

% 6VesNCycles
[ang_locs_tagmode,final_encoding_matrix] = hw_fun_IterEncoding_FixedHadamard(VesLocs3D,AngTowardsCor,AngTowardsSag,off_res,...
        number_of_iterations,modmat_phase,modmat_lookup,opts); clc;

fprintf('Please enter the SubID (an integer): ');
while true
    SubID = input('');
    if isnan(SubID) || mod(SubID, 1) ~= 0
        fprintf('The input is not an integer. Please try again.\n');
    else
        break;
    end
end

fprintf('Please enter the SessID (an integer): ');

while true
    SessID = input('');
    if isnan(SessID) || mod(SessID, 1) ~= 0
        fprintf('The input is not an integer. Please try again.\n');
    else
        break;
    end
end

if if_fft_based
    save_path = [workingDir,filesep,'IOES_6Ves',num2str(Hadamard_order),'Cycles','_',vendor,'_isUnipolar',num2str(isUnipolar),'_Sub',num2str(SubID),'_Sess',num2str(SessID)];
else
    save_path = [workingDir,filesep,'MOGEN_6Ves',num2str(Hadamard_order),'Cycles','_',vendor,'_isUnipolar',num2str(isUnipolar),'_Sub',num2str(SubID),'_Sess',num2str(SessID)];
end
if exist(save_path,'dir')==7
    rmdir(save_path,'s');
end
mkdir(save_path);

writematrix(ang_locs_tagmode,[save_path,filesep,'VEPCASL_OES_Setup.txt']);
writematrix(VesLocs3D,[save_path,filesep,'VesLocs3D.txt']); % order-RICA,RVA,LICA,LVA,RECA,LECA

rmpath(util_path);