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

function [ cycles, duration ] = find_cycles( mag, est_period)
%find_cycles - Split signal into gait cycles
% mag        - Acceleration magnitude
% est_period - Estimated period of gait cycles
% cycles     - n x 2 matrix of (start, end) indices of gait cycles
% duration   - Number of samples in cycle

% Find peaks in mag
ld = find([0; mag(2:end)-mag(1:end-1)] > 0);
rd = find([mag(1:end-1)-mag(2:end); 0] > 0);
peaks = intersect(ld, rd);

% Kill lower peaks < 20% est_period from another peak
kill = zeros(size(peaks));
for i=1:numel(peaks)
    lo = max(floor(peaks(i) - est_period / 5), 1);
    hi = min(ceil(peaks(i) + est_period / 5), numel(mag));
    if (max(mag(lo:hi)) > mag(peaks(i)))
        kill(i) = 1;
    end
end
peaks = peaks(kill == 0);

% Keep only peaks with other peaks in est_period * (1 +- 20%)
kill = zeros(size(peaks));
for i=1:numel(peaks)
    t1 = abs(peaks-peaks(i));
    t2 = abs(t1 - est_period);
    t3 = min(t2);
    if (t3 > (est_period * 20 / 100))
        kill(i) = 1;
    end
end
peaks = peaks(kill == 0);

cycles = zeros(numel(peaks), 2, 'int32');
duration = zeros(numel(peaks), 1);
ci = 0;
for i=1:numel(peaks)
    cycles(ci+1,1) = peaks(i);
    % Find another paek est_period * (1 += 20%)
    lo = peaks(i) + floor(est_period * 0.8);
    hi = peaks(i) + ceil(est_period * 1.2);
    pidx = find((peaks >= lo) .* (peaks <= hi));
    if (numel(pidx) == 1)
        cycles(ci+1, 2) = peaks(pidx);
        ci = ci + 1;
        duration(ci) = cycles(ci, 2) - cycles(ci, 1);
    end
end
cycles = cycles(1:ci, :);
duration = duration(1:ci);
end
