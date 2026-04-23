function [ minValue,closestIndex ] = findVectorClosestValue( vTarget,vTimePoints )

% Retrieve closest values to a predefined vector and its position within
% the vector.
% 
% vTarget: vector that contains the set of values to find
%
% vTimePoints: set of values that may not match elements of vTarget 
%

%%

minValue = [];
closestIndex = [];

for ind=1:length(vTimePoints)
    [ ~,pos ] = min( abs(vTarget - vTimePoints(ind)) );

    minValue = [ minValue; vTarget(pos) ];
    closestIndex = [ closestIndex; pos ];
end

end