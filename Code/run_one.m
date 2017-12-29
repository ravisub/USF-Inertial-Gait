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

function run_one (datadir, probelen, use_linacc, use_gyro, rotate)
% run_one    - Runs experiment specified by the data direectory
% data_dir   - Path to directory containing sensor csv files
% probelen   - Duration of probe in seconds
% use_linacc - Boolean, If Linear acceleration should be used for matching
% use_gyro   - Boolean, If Gyroscope should be used for matching
% rotate     - Boolean, If rotation should be performed before matching

%% Load csv files into gallery and probes
tic
[gallery, probes] = load_gallery_probes([datadir '\'], probelen);
disp(['Loaded gallery and probes in ' num2str(toc) ' seconds'])

%% Score the probes against gallery
tic
score = gait_score(gallery, probes, use_linacc, use_gyro, rotate);
disp(['Scored in ' num2str(toc) ' seconds'])

%% Create Normalized score each probe against gallery
tic
score(:,5) = -score(:,5);
scoren = score;
up = unique(scoren(:,2));
for i = 1:numel(up)
    ms = mean(scoren(scoren(:,2)==up(i),5));
    ss = std(scoren(scoren(:,2)==up(i),5));
    scoren(scoren(:,2)==up(i), 5) = (scoren(scoren(:,2)==up(i), 5) - ms) / ss;
end
clear i ms ss up
disp(['Nornalized scores in ' num2str(toc) ' seconds'])

%% Compute ROC & confidence intervals (Unnormalized)
tic
rocin(:,2) = score(:,5);
rocin(score(:,2) == score(:,4), 1) = 1;
rocin(score(:,2) ~= score(:,4), 1) = 0;
dlmwrite('rocin.txt', rocin, ' ');
system('.\CVRLROC.exe -s rocin.txt -r rocout.m -n 100 -b -l nn');
rocout
[a b] = min(abs(meanFA_nn + meanTA_nn - 1));
eernn = meanFA_nn(b);
[a b] = min(abs(CIFA_low_nn + CITA_high_nn - 1));
eernn_low = CIFA_low_nn(b);
[a b] = min(abs(CIFA_high_nn + CITA_low_nn - 1));
eernn_high = CIFA_high_nn(b);
clear a b
disp(['Computed unnornalized ROC in ' num2str(toc) ' seconds'])

%% Compute ROC & confidence intervals (Normalized)
tic
rocin(:,2) = scoren(:,5);
rocin(scoren(:,2) == scoren(:,4), 1) = 1;
rocin(scoren(:,2) ~= scoren(:,4), 1) = 0;
dlmwrite('rocin.txt', rocin, ' ');
system('.\CVRLROC.exe -s rocin.txt -r rocout.m -n 100 -b -l yn');
rocout
[a b] = min(abs(meanFA_yn + meanTA_yn - 1));
eeryn = meanFA_yn(b);
[a b] = min(abs(CIFA_low_yn + CITA_high_yn - 1));
eeryn_low = CIFA_low_yn(b);
[a b] = min(abs(CIFA_high_yn + CITA_low_yn - 1));
eeryn_high = CIFA_high_yn(b);
clear a b
disp(['Computed nornalized ROC in ' num2str(toc) ' seconds'])

%% Save results in .mat file
tic
[~, expt] = fileparts(datadir);
save(['..\Results\' expt '_p' num2str(probelen, '%02d'), '_linacc_', num2str(use_linacc), '_gyro_', num2str(use_gyro), '_rotate_' num2str(rotate) '.mat']);
disp(['Saved in ' num2str(toc) ' seconds'])
