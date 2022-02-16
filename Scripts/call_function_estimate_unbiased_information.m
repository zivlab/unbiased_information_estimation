%% Main script for correcting the upward bias in the naive calculation of information content:

% Here, we use two methods to correct the upward bias in the naive calculation of information content,
% seeking a bias-free estimation of the information carried by neuronal
% activity about a given encoded variable:
% -----------------------------------------

% 1. The scaled shuffle reduction (SSR) method - estimating the ratio between the
% bias in the naive and shuffle information and subtracting from the naïve information
% the shuffle information scaled by the estimated bias ratio.

% 2. The bounded asymptotic extrapolation (BAE) method - fitting the function of how the information changes with sample size and
% extrapolating it to infinity.

% *The naive estimation of information content is based either on the Skaggs information index (SI) or on the mutual information (MI).
% **The bias correction methods can be modified to accomodate additional information theoretic measures.
% ***The obtained results are compared against the shuffle reduction (SR) and unbounded asymptotic extrapolation (AE) methods.

% The inputs for this script should include:
% -------------------------------------------
% 1. spike_train - a matrix of size TxN, where T is the number of time bins
% and N is the number of neurons. Each element is the spike count of a
% given neuron in a given time bin.
% 2. stimulus_trace - a vector of size T, where each element is the stimulus value in a
% given time bin.

% *When applied to quantify the spatial tuning of place cells, the common practice is to first filter out 
% time bins durning which the animal is not running (i.e., analyze only bins with speed > threshold) 

% Outputs:
%----------
% General settings (dt, number of shuffles, subsampling repetitions, etc.) 
% Firing statistics (average rates, average active time bins, fraction of significantly tuned cells, etc.)
% Estimated information (SI and MI based on the SSR and BAE methods)

clc
close all

%% Defining all the settings and parameters for calling the unbiased_information_estimation function:

settings=struct;
% Choosing which information theoretic measures to estimate:
estimate_SI_bit_spike=1; % Choose to estimate SI (bit/spike) by setting the value to 1 (0 otherwise)
estimate_SI_bit_sec=1; % Choose to estimate SI (bit/sec) by setting the value to 1 (0 otherwise)
estimate_MI=0; % Choose to estimate MI by setting the value to 1 (0 otherwise)
measures_to_estimate=[estimate_SI_bit_spike,estimate_SI_bit_sec,estimate_MI];
settings.measures_to_estimate=measures_to_estimate; 

% General parameters:
settings.dt=0.05; % time bin in units of seconds

% Focusing only on sufficiently active cells:
settings.active_bins_threshold=10; % minimal number of time bins in which the cell is defined active - less than 10 active time bins leads to an inaccurate estimation
settings.firing_rate_threshold=0; % in spike/sec.  Default value is 0, but you can choose to add an average firing rate threshold

% Settings for the tuning signficance test:
settings.estimate_only_significantly_tuned_cells=1; % 1 if estimation is performed only for significanty tuned cells (0 otherwise)
settings.shuffle_type='cyclic'; % permutations used for the significance test can be either 'cyclic' or 'random'
settings.num_shuffles=1000;
settings.tuning_significance_threshold=0.05;

% Setting for the subsampling procedure:
settings.subsampling_repetitions=500; % number of repetitions in the subsampling of the data
settings.subsample_fraction=(0.1:0.1:1); % subsamples with size of different fractions of the data

settings.plot_results=0; % 1 for plotting the results (0 otherwise)
settings.save_figures=0; % 1 for saving the figures (0 otherwise)
settings.figures_directory=[]; % if plotting results then set the path for saving the figures

%% loading the data and calling the unbiased_information_estimation function: 

data_path='D:\dev\Bias correction\unbiased_information_estimation\Real data\Sample data\';
temp_spike_train=load(fullfile(data_path,'spike_train.mat'));
spike_train=temp_spike_train.spike_train;
temp_stimulus_trace=load(fullfile(data_path,'stimulus_trace.mat'));
stimulus_trace=temp_stimulus_trace.stimulus_trace;

% performing the estimation of information content:
unbiased_information_estimation_results=estimate_unbiased_information(spike_train,stimulus_trace,settings);
