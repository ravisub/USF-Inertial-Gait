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

%%
% This script is used to run many experiments in batch mode by calling the
% run_one script for each experiment.
% Experiment Same_day_folder Days_apart_folder
% 1          probes          repeats
% 2          probes_13-15    repeats_13-15
% 3          probes_16-18    repeats_16-18
% 5          probes_13       repeats_13
% 6          probes_16       repeats_6

%% 
run_one('..\Data\probes', 5, 1, 1, 0);
run_one('..\Data\probes', 5, 1, 1, 1);
run_one('..\Data\probes_13-15', 5, 1, 1, 0);
run_one('..\Data\probes_13-15', 5, 1, 1, 1);
run_one('..\Data\probes_16-18', 5, 1, 1, 0);
run_one('..\Data\probes_16-18', 5, 1, 1, 1);
run_one('..\Data\probes_13', 5, 1, 1, 0);
run_one('..\Data\probes_13', 5, 1, 1, 1);
run_one('..\Data\probes_16', 5, 1, 1, 0);
run_one('..\Data\probes_16', 5, 1, 1, 1);
% 
%%
run_one('..\Data\repeats', 5, 1, 1, 0);
run_one('..\Data\repeats', 5, 1, 1, 1);
run_one('..\Data\repeats_13-15', 5, 1, 1, 0);
run_one('..\Data\repeats_13-15', 5, 1, 1, 1);
run_one('..\Data\repeats_16-18', 5, 1, 1, 0);
run_one('..\Data\repeats_16-18', 5, 1, 1, 1);
run_one('..\Data\repeats_13', 5, 1, 1, 0);
run_one('..\Data\repeats_13', 5, 1, 1, 1);
run_one('..\Data\repeats_13', 5, 1, 1, 1);
run_one('..\Data\repeats_16', 5, 1, 1, 0);
run_one('..\Data\repeats_16', 5, 1, 1, 1);
