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

function [ gallery, probes ] = load_gallery_probes( datadir, probelen )
% load_gallery_probes - Loads data into gallery and probes
% data_dir   - Path to directory containing sensor csv files
% probelen   - Duration of probe in seconds
%
% data_dir should contain Walks.txt that specifies how to split each csv
% file into gallery and probes.

probelen = probelen * 100; % Convert to # of samples

l = dir([datadir '*.csv']);
f = fopen([datadir 'Walks.txt']);
gallery_idx = 0;
probe_idx = 0;
tgp_struct = struct('id', 0, 'activity', 0, 'range', {}, 'us_sig', {});
tgallery = tgp_struct;
tprobes = tgp_struct;
for i = 1:numel(l)
    id = str2num(l(i).name(2:4));
    activity = str2num(l(i).name(6:7));
    fid = fscanf(f, '%d', 1);
    fact = fscanf(f, '%d', 1);
    fnsub = fscanf(f, '%d', 1);
    if ((str2num(l(i).name(1:4)) ~= fid) | (activity ~= fact))
        disp(['Mismatch Line - ' num2str(i) ' id - ' num2str(id) ' fid - ' num2str(fid) ' activity - ' num2str(activity) ' fact - ' num2str(fact)])
    end
    data = csvread([datadir l(i).name], 1, 0);
    for j = 1:fnsub
        range = fscanf(f, '%d', 2);
        seq = tgp_struct;
        seq(1).id = id;
        seq(1).activity = activity;
        seq(1).range = range;
        raw_sig.gyro = data(range(1):range(2), [12 9:11]);
        raw_sig.linacc = data(range(1):range(2), [16 13:15]);
        seq(1).us_sig = uniform_sample(raw_sig, 0.01);
        gorp = 'p'; % Initialize to probe
        if (strfind(datadir, 'repeats'))
            if (l(i).name(1) == '1')
                if (j == 1)
                    gorp = 'g';
                else
                    continue;
                end
            end
        else
            if (((activity == 13) | (activity == 16)) & (j == 1))
                gorp = 'g';
            end
        end
        if (gorp == 'g')
            gallery_idx = gallery_idx + 1;
            tgallery(gallery_idx) = seq(1);
        else
            probe_idx = probe_idx + 1;
            tprobes(probe_idx) = seq(1);
        end
    end
end
fclose(f);
clear activity data f fact fi fid fnsub gallery_idx i id j l probe_idx range raw_sig seq

%% Process gallery
for i=1:numel(tgallery)
    period = acperiod(tgallery(i).us_sig.linacc);
    tgallery(i).us_sig.mag = sqrt(sum(tgallery(i).us_sig.linacc .^ 2, 2));
    cycle_marks = find_cycles(tgallery(i).us_sig.mag, period);
    cycles = [tgallery(i).us_sig.mag tgallery(i).us_sig.linacc tgallery(i).us_sig.gyro];
    resampled_cycles = resample_cycles(cycles, cycle_marks, 100);
    tgallery(i).cycle_marks = cycle_marks;
    tgallery(i).cycles = cycles;
    tgallery(i).resampled_cycles = resampled_cycles;
end

% Combine Pocket and Holster
gp_struct = struct('id', 0, 'activity', 0, 'cycle_marks', {}, 'cycles', {}, 'resampled_cycles', {});
gallery = gp_struct;
uid = unique([tgallery.id]);
for i=1:numel(uid)
    ttt = find([tgallery.id] == uid(i));
    gallery(i).id = tgallery(ttt(1)).id;
    gallery(i).activity = tgallery(ttt(1)).activity;
    maxt = 0;
    gallery(i).cycles = tgallery(ttt(1)).cycles;
    gallery(i).resampled_cycles = tgallery(ttt(1)).resampled_cycles;
    gallery(i).cycle_marks = tgallery(ttt(1)).cycle_marks + maxt;
    maxt = maxt + size(gallery(i).cycles, 1);
    for j=2:numel(ttt)
        gallery(i).cycles = [gallery(i).cycles; tgallery(ttt(j)).cycles];
        gallery(i).resampled_cycles = [gallery(i).resampled_cycles; tgallery(ttt(j)).resampled_cycles];
        gallery(i).cycle_marks = [gallery(i).cycle_marks; tgallery(ttt(j)).cycle_marks + maxt];
        maxt = maxt + size(gallery(i).cycles, 1);
    end
end
clear cycle_marks cycles i j maxt period resampled_cycles ttt uid

%% Process probes
probes = gp_struct;
probe_idx = 0;
for i=1:numel(tprobes)
    tprobes(i).us_sig.mag = sqrt(sum(tprobes(i).us_sig.linacc .^ 2, 2));
    lim = floor(numel(tprobes(i).us_sig.mag) / probelen) * probelen;
    for j=1:probelen:lim
        period = acperiod(tprobes(i).us_sig.linacc(j:j+probelen-1,:));
        if (numel(period) > 0)
            probe_idx = probe_idx + 1;
            probes(probe_idx).tp = i;
            probes(probe_idx).range = [j j+probelen-1];
            probes(probe_idx).id = tprobes(i).id;
            probes(probe_idx).activity = tprobes(i).activity;
            probes(probe_idx).us_sig.mag = tprobes(i).us_sig.mag(j:j+probelen-1,:);
            probes(probe_idx).us_sig.linacc = tprobes(i).us_sig.linacc(j:j+probelen-1,:);
            probes(probe_idx).us_sig.gyro = tprobes(i).us_sig.gyro(j:j+probelen-1,:);
            cycle_marks = find_cycles(tprobes(i).us_sig.mag(j:j+probelen-1), period);
            cycles = [tprobes(i).us_sig.mag(j:j+probelen-1,:) tprobes(i).us_sig.linacc(j:j+probelen-1,:) tprobes(i).us_sig.gyro(j:j+probelen-1,:)];
            resampled_cycles = resample_cycles(cycles, cycle_marks, 100);
            probes(probe_idx).cycle_marks = cycle_marks;
            probes(probe_idx).cycles = cycles;
            probes(probe_idx).resampled_cycles = resampled_cycles;
        end
    end
    % Take care of remaining length < probelen
    period = acperiod(tprobes(i).us_sig.linacc(lim+1:end,:));
    if (numel(period) > 0)
        probe_idx = probe_idx + 1;
        probes(probe_idx).tp = i;
        probes(probe_idx).range = [lim+1 numel(tprobes(i).us_sig.mag)];
        probes(probe_idx).id = tprobes(i).id;
        probes(probe_idx).activity = tprobes(i).activity;
        probes(probe_idx).us_sig.mag = tprobes(i).us_sig.mag(lim+1:end,:);
        probes(probe_idx).us_sig.linacc = tprobes(i).us_sig.linacc(lim+1:end,:);
        probes(probe_idx).us_sig.gyro = tprobes(i).us_sig.gyro(lim+1:end,:);
        cycle_marks = find_cycles(tprobes(i).us_sig.mag(lim+1:end), period);
        cycles = [tprobes(i).us_sig.mag(lim+1:end,:) tprobes(i).us_sig.linacc(lim+1:end,:) tprobes(i).us_sig.gyro(lim+1:end,:)];
        resampled_cycles = resample_cycles(cycles, cycle_marks, 100);
        probes(probe_idx).cycle_marks = cycle_marks;
        probes(probe_idx).cycles = cycles;
        probes(probe_idx).resampled_cycles = resampled_cycles;
    end
end

clear cycle_marks cycles i j lim period probe_idx resampled_cycles
end
