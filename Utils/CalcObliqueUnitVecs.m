% Calculates unit vectors in the "AP" and "LR" directions for oblique
% labelling slabs with angles defined as per the Siemens UI
%
% [LRUnitVec APUnitVec] = CalcObliqueUnitVecs(AngTowardsCor,AngTowardsSag)

function [LRUnitVec, APUnitVec] = CalcObliqueUnitVecs(AngTowardsCor,AngTowardsSag)

% Start with standard Unit vectors
LRUnitVec = [1,0,0]';
APUnitVec = [0,1,0]';

% Calculate rotation towards coronal
RotMatTowardsCor = RotMat3D_TO([1 0 0],todeg2rad(-AngTowardsCor));

% Calculate the new AP unit vector
APUnitVec = RotMatTowardsCor * APUnitVec;

% The next rotation towards sagittal doesn't affect the new AP unit vector
% so must be about the AP unit vector itself
RotMatTowardsSag = RotMat3D_TO(APUnitVec,todeg2rad(-AngTowardsSag));

% Calculate the new LR unit vector
LRUnitVec = RotMatTowardsSag * LRUnitVec;
