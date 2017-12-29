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

function [ score ] = compare( pcycles, gcycles, method )
%compare Compares probe and gallery using method
% pcycles - probe signal
% gcycles - gallery cycle
% method - one of 'L1','L2', 'Tanimoto', 'Dot', 'NCC'
if (strcmpi(method, 'L1'))
    score = sum(abs(pcycles(:)-gcycles(:)));
    return
end
if (strcmpi(method, 'L2'))
    score = sum((pcycles(:)-gcycles(:)).^2);
    return
end
if (strcmpi(method, 'Tanimoto'))
    ab = sum(pcycles(:) .* gcycles(:));
    aa = sum(pcycles(:) .^ 2);
    bb = sum(gcycles(:) .^ 2);
    score = -ab ./ (aa + bb - ab);
    return
end
if (strcmpi(method, 'Dot'))
    ab = sum(pcycles(:) .* gcycles(:));
    aa = sum(pcycles(:) .^ 2);
    bb = sum(gcycles(:) .^ 2);
    score = -ab / sqrt(aa * bb);
    return
end
if (strcmpi(method, 'NCC'))
    pcycles = pcycles - repmat(mean(pcycles), size(pcycles,1), 1);
    pcycles = pcycles(:) / sqrt(sum(pcycles(:) .* pcycles(:)));
    gcycles = gcycles - repmat(mean(gcycles), size(gcycles,1), 1);
    gcycles = gcycles(:) / sqrt(sum(gcycles(:) .* gcycles(:)));
    score = -sum(pcycles .* gcycles);
    return
end
error(['Unknown comparison method (' method ') - must be one of L1 L2 Tanimoto NCC or DTW'])
end
