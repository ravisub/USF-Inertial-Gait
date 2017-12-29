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

function [ us_signals ] = uniform_sample( raw_signals, sp )
% uniform_sample - Lines up and resamples the raw_signals at uniform rate
% raw_signals - Linear acceleration and Gyroscope data with timestamps
% sp          - Time between samples 0.01 for 100 Hz

%% Find common domain for signals
min_ts = max([min(raw_signals.gyro(:,1)) min(raw_signals.linacc(:,1))]);
max_ts = min([max(raw_signals.gyro(:,1)) max(raw_signals.linacc(:,1))]);
% Set lowest timestamp to zero and convert to seconds
raw_signals.gyro(:,1) = (raw_signals.gyro(:,1) - min_ts) / 1e9;
raw_signals.linacc(:,1) = (raw_signals.linacc(:,1) - min_ts) / 1e9;
max_ts = (max_ts - min_ts) / 1e9;
min_ts = (min_ts - min_ts) / 1e9; % Set to zero

%% Resample
rs_ts = [min_ts:sp:(max_ts-sp)];
us_signals.ts = rs_ts;
us_signals.gyro = uniform_sample_one(raw_signals.gyro, rs_ts);
us_signals.linacc = uniform_sample_one(raw_signals.linacc, rs_ts);
us_signals.duration = rs_ts(end);
end

function [resampled] = uniform_sample_one (raw_signal, rs_ts)
i = 1;
j = 1;
for t = rs_ts
    while(raw_signal(i,1) <= t)
        i = i + 1;
    end
    a = t - raw_signal(i-1,1);
    b = raw_signal(i, 1) - t;
    resampled(j, :) = (b * raw_signal(i-1, 2:end) + a * raw_signal(i, 2:end)) / (a + b);
    j = j + 1;
end
end
