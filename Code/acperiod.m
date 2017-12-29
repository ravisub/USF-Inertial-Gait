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

function [ period ] = acperiod( sig )
% acperiod - Estimates period using auto-correlation
nominal_period = 100; % 1 second
len = size(sig, 1);
if (len == 0)
    period = [];
    return
end
for i=1:ceil(len / 2)
    sh = i-1;
    ssig = sig([sh+1:end 1:sh],:);
    dp = sig .* ssig;
    ac(i) = sum(dp(:));
end
% Finds peaks in ac
ld = find([0 ac(2:end)-ac(1:end-1)] > 0);
rd = find([ac(1:end-1)-ac(2:end) 0] > 0);
peaks = intersect(ld, rd);
% Kill lower peaks < 20% est_period from another peak
kill = zeros(size(peaks));
for i=1:numel(peaks)
    lo = max(floor(peaks(i) - nominal_period / 5), 1);
    hi = min(ceil(peaks(i) + nominal_period / 5), numel(ac));
    if (max(ac(lo:hi)) > ac(peaks(i)))
        kill(i) = 1;
    end
end
peaks = peaks(kill == 0);
% Pick peak closest to nominal_period
[~, np] = min(abs(peaks-nominal_period));
period = peaks(np);
end
