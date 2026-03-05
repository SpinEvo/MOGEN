% Calculates 2D vessel locations from 3D vessel locations
function VesLocs2D = Convert3DTo2DVesLocs(VesLocs3D, LRUnitVec, APUnitVec)

VesLocs3D = VesLocs3D'; %convert column position vectors in to row vectors
VesLocs2D = zeros(size(VesLocs3D,1),2);
% Loop through and calculate the amount each vessel lies in the direction
% of the two unit vectors via dot products
for ii=1:size(VesLocs3D,1)
   VesLocs2D(ii,:) = [sum(VesLocs3D(ii,:).*LRUnitVec(:)') sum(VesLocs3D(ii,:).*APUnitVec(:)')];
end

VesLocs2D = VesLocs2D';