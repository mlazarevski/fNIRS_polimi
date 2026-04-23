function [ T,pMovingNew ] = computeAffineTransformFromPointcloud( pMoving,pFixed )

% This function computes the affine transform T with respect to two cloud
% of corresponding points. 
% 
% INPUT:
% - pMoving: cloud of points - moving points
% - pFixed: cloud of points - fixed points
% 
% OUTPUT:
% - T: affine transform
% - pMovingNew: pMoving after the application of T transform
%

%%

pMoving = [ pMoving ones(size(pMoving, 1),1) ]';
pFixed = [ pFixed ones(size(pFixed, 1),1) ]';
T = pFixed / pMoving;

pMoving = pMoving';
pFixed = pFixed';
pMoving = pMoving( :,1:3 );
pFixed = pFixed( :,1:3 );
pMovingNew = xform_apply(pMoving, T);

end

