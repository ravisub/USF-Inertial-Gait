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

function [ score ] = gait_score( gallery, probes, use_linacc, use_gyro, rotate )
% gait_score - Computes score for each probe-gallery pair
% gallery    - Gallery
% probes     - Probes
% use_linacc - Use Linear Acceleration to compute scoree
% use_gyro   - Use Gyroscope to compute scores
% rotate     - Rotate probe before computing score

si = 1;
score = zeros(numel(probes) * numel(gallery), 5);
% atscores = zeros(0, 5);
for p=1:numel(probes)
    for g=1:numel(gallery)
        score(si, 1) = p;
        score(si, 2) = probes(p).id;
        score(si, 3) = g;
        score(si, 4) = gallery(g).id;
        tscores = ones(size(probes(p).resampled_cycles,1) * size(gallery(g).resampled_cycles,1), 1) * Inf;
        if (numel(tscores) == 0)
            continue;
        end
        ti = 1;
        for pc = 1:size(probes(p).resampled_cycles,1)
            pcycles = zeros(size(probes(p).resampled_cycles,2),0);
            if (use_linacc)
                pcycles(:,end+1:end+3) = squeeze(probes(p).resampled_cycles(pc,:,2:4));
            end
            if (use_gyro)
                pcycles(:,end+1:end+3) = squeeze(probes(p).resampled_cycles(pc,:,5:7));
            end
            for gc = 1:size(gallery(g).resampled_cycles,1)
                gcycles = zeros(size(gallery(g).resampled_cycles,2),0);
                if (use_linacc)
                    gcycles(:,end+1:end+3) = squeeze(gallery(g).resampled_cycles(gc,:,2:4));
                end
                if (use_gyro)
                    gcycles(:,end+1:end+3) = squeeze(gallery(g).resampled_cycles(gc,:,5:7));
                end
                if (rotate)
%                    rpcycles = rotate_distance(squeeze(probes(p).resampled_cycles(pc,:,2:4)), squeeze(gallery(g).resampled_cycles(gc,:,2:4)), pcycles); % Use linacc to rotate
                    rpcycles = rotate_distance(pcycles, gcycles, pcycles); % Use all signal to rotate
                else
                    rpcycles = pcycles;
                end
                tscores(ti) = compare(rpcycles, gcycles, 'Tanimoto');
%                atscores(end+1,:) = [p pc g gc tscores(ti)];
                ti = ti + 1;
            end
        end
        score(si,5) = min(tscores);
        si = si + 1;
    end
    if (mod(p,100)==0)
        ['Scored ' num2str(p) ' of ' num2str(numel(probes))]
    end
end
% save(['atscores_' num2str(rotate) '.mat'], atscores);
score = score(1:si-1,:);
clear si p g tscores ti pc pcycles gc gcycles
end
