%% generate data

%Please be aware that generating each data file is time/CPU intensive. 
%Sections (below) are ordered approximately by time to run ~per file~, in 
%ascending order.
%We recommend generating only the data required by evaluating the appropriate section.
%Reducing ntrials and nsamples by a factor of ten will give adequate
%interim results.

%Additionally, please note that the split data for fig5(d-f) can be found in the
%data.bris respository linked to the NatComm22 paper, in subfolder '3N curve SAT param 
%split'. 
%We do not provide code for splitting the data and converting to .txt here.
%The data is taken from 3N_curve_data_dense.mat, should you wish to do this
%yourself.

%% add utils to the path

utils_folder = fullfile(cd, '..', 'utils');

if exist(utils_folder, 'dir')
    addpath(utils_folder);
else
    error('Please download the utils folder and contents.');
end

%% curve parameterisation by N

folder = fullfile(cd, '..', 'data', 'N-choice curve');
if ~exist(folder, 'dir')
    mkdir(folder)
end

nt = 100; 
ntrials = 1000000;
nsamples = 100000;

% curve parameterisation by N: required for fig9a, fig9b
generate_curve_data_cell(3, fullfile('N-choice curve', '3N_curve_data.mat'), nt, ntrials, nsamples);
generate_curve_data_cell(5, fullfile('N-choice curve', '5N_curve_data.mat'), nt, ntrials, nsamples);
generate_curve_data_cell(7, fullfile('N-choice curve', '7N_curve_data.mat'), nt, ntrials, nsamples);
generate_curve_data_cell(9, fullfile('N-choice curve', '9N_curve_data.mat'), nt, ntrials, nsamples);

%% curve band analysis

folder = fullfile(cd, '..', 'data', 'band analysis');
if ~exist(folder, 'dir')
    mkdir(folder)
end

ntrials = 1000000;
nsamples = 100000;

% curve band analysis: required for fig6
generate_banded_curve_data(3, fullfile('band analysis', '3N_curve_data_banded.mat'), ntrials, nsamples)

%% cuve, oscil, and power parameterisations for 3N

folder = fullfile(cd, '..', 'data', '3N parameterisations');
if ~exist(folder, 'dir')
    mkdir(folder)
end

ntrials = 1000000;
nsamples = 100000;
nt = 50;

% oscil parameterisation dense: required for fig4b, fig8cfi
generate_oscil_data(3, fullfile('3N parameterisations', '3N_oscil_data.mat'), nt, ntrials, nsamples);

% power parameterisation dense: required for fig4c, fig8beh
generate_power_data(3, fullfile('3N parameterisations', '3N_power_data.mat'), nt, ntrials, nsamples);

nt = 500;
% curve parameterisation dense: required for fig4a, fig5(all), fig8adg
generate_curve_data(3, fullfile('3N parameterisations', '3N_curve_data_dense.mat'), nt, ntrials, nsamples);

