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

function [ r_cycles ] = resample_cycles( raw, cycles, N )
% resample_cycles   Create N samples
% raw            - Raw acceleration data vertical, lateral, forward
% cycles         - n x 3 matrix of (start, end, status) indices of gait cycles
%                - status 0, 1  - bad, good
% N              - number of resamples per cycle
% r_cycles       - Resampled cycle

r_cycles = zeros(size(cycles,1), N+2, size(raw,2));
for i=1:size(cycles,1)
    raw_cycle = raw(cycles(i,1):cycles(i,2), :);
    r_cycle = zeros(N+2, size(raw,2));
    % Normalize timestamps
    lll = 1;
    hhh = size(raw_cycle,1);
    n_ts = ([lll:hhh] - lll) * (N+1) / (hhh - lll);

    % Time normalization & Resampling
    right = 1;
    for j=0:N
        while (n_ts(right) <= j) % Fine right point
            right = right + 1;
        end
        left = right - 1;
        r_cycle(j+1,:) = raw_cycle(left,:) + (raw_cycle(right,:) - raw_cycle(left,:)) * (j - n_ts(left)) / (n_ts(right) - n_ts(left));
    end
    r_cycle(N+2,:) = raw_cycle(end,:);
    r_cycles(i, :, :) = r_cycle;
end
end
