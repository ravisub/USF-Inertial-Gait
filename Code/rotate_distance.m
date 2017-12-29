%%
% Copyright (c) 2018 by University of South Florida, Tampa, FL, USA.
% Use, or copying without permission prohibited.
% PERMISSION TO USE
% In transmitting this software, permission to use for research and
% educational purposes is hereby granted. This software may be copied for
% archival and backup purposes only. This software may not be transmitted
% to a third party without prior permission of the copyright holder.
% This permission may be granted only by Prof. Sudeep Sarkar of University
% of South Florida (sarkar@mail.usf.edu).
% Acknowledgment as appropriate is respectfully requested.

%%
% Author
% Ravichandran Subramanian
% Department of Computer Science and Engineering
% University of South Florida

function [ rprobes] = rotate_distance( pref, gref, probes )
% rotate_distance - Rotates probe to match gallery using Kabsch algorithm
% pref            - Probe signal to match gallery
% gref            - Gallery signal to be matched
% probes          - Probe signals to be rotated
pref3 = zeros(0,3);
gref3 = zeros(0,3);
for i=1:3:size(pref,2)
    pref3 = [pref3; pref(:,i:i+2)];
    gref3 = [gref3; gref(:,i:i+2)];
end
a = pref3' * gref3;
[v s w] = svd(a);
d = sign(det(w * v'));
u = w * ([1 0 0; 0 1 0; 0 0 d]) * v';
for i=1:3:size(probes,2)
    rprobes(:,i:i+2) = (u * probes(:,i:i+2)')';
end
end
