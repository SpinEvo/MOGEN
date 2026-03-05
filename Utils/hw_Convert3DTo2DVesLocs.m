function VesLocs2D = hw_Convert3DTo2DVesLocs(VesLocs3D, LRUnitVec, APUnitVec)

% Loop through and calculate the amount each vessel lies in the direction
% of the two unit vectors via dot products
for ii=1:size(VesLocs3D,1)
   VesLocs2D(ii,:) = [sum(VesLocs3D(ii,:).*LRUnitVec(:)') sum(VesLocs3D(ii,:).*APUnitVec(:)')];
end