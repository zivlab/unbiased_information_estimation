%% Main script for correcting the upward bias in the naive calculation of information content:

% Here, we use two methods to correct the upward bias in the naive calculation of information content,
% seeking a bias-free estimation of the information carried by neuronal
% activity about a given encoded variable:
% -----------------------------------------

% 1. The scaled shuffle reduction (SSR) method - estimating the ratio between the
% bias in the naive and shuffle information and subtracting from the naïve information
% the shuffle information scaled by the estimated bias ratio.

% 2. The bounded extrapolation (BAE) method - fitting the function of how the information changes with sample size and
% extrapolating it to infinity.

% *The naive estimation of information content is based either on the Skaggs information index (SI) or on the mutual information (MI).
% **The bias correction methods can be modified to accomodate additional information theoretic measures.
% ***The obtained results are compared against the shuffle reduction (SR) and unbounded extrapolation (AE) methods.

% The inputs for this script should include:
% -------------------------------------------
% 1. spike_train - a matrix of size TxN, where T is the number of time bins
% and N is the number of neurons. Each element is the spike count of a
% given neuron in a given time bin.
% 2. stimulus_trace - a vector of size T, where each element is the stimulus value in a
% given time bin.

% Outputs:
%----------
% General parameters (dt, number of shuffles, subsampling repetitions, etc.) and
% data statistics (average rates, average active time bins, fraction of significant cells, etc.)
% Estimated information (SI and MI based on the SSR and BAE methods)

clc
close all

%% Choose the information measures you would like to estimate by setting their values to 1 (0 otherwise):

estimate_SI_bit_spike=1;
estimate_SI_bit_sec=1;
estimate_MI=1;
measures_to_estimate=[estimate_SI_bit_spike,estimate_SI_bit_sec,estimate_MI];

%% Step 1 - load the spike train and stimulus trace:

data_path='D:\dev\Bias correction\unbiased_information_estimation\Simulated data\Sample data\';
temp_spike_train=load(fullfile(data_path,'simulated_spike_train.mat'));
simulated_spike_train=temp_spike_train.simulated_spike_train;
temp_stimulus_trace=load(fullfile(data_path,'stimulus_trace.mat'));
stimulus_trace=temp_stimulus_trace.stimulus_trace;
temp_true_information=load(fullfile(data_path,'true_information.mat'));
true_information=temp_true_information.true_information;
mkdir(data_path,'Sample data - results')
mkdir(fullfile(data_path,'Sample data - results'),'Figures')
figures_directory=[data_path '\Sample data - results\Figures'];

%% Step 2: Identifying significantly modulated cells by comparing the naive information with shuffles:

% Setting the parameters and which cells to analyze:
dt=1/20; % in units of seconds
shuffle_type='cyclic'; % permutations can be either 'cyclic' or 'random'
active_bins_threshold=10; % minimal number of time bins in which the cell is defined active - less than 10 active time bins leads to an inaccurate estimation
firing_rate_threshold=0; % in spike/sec.  Default value is 0, but you can choose to add an average firng rate threshold
num_shuffles=1000;
significance_threshold=1; % for determining which cells are significantly tuned to the encoded variable

% Computing the rate maps for shuffled spike trains:
[rate_maps,average_firing_rates,normalized_states_distribution]=compute_rate_maps(simulated_spike_train,stimulus_trace,dt);
active_bins=sum(simulated_spike_train>0);
active_cells=find(active_bins>=active_bins_threshold & average_firing_rates>firing_rate_threshold);
shuffled_spike_trains=shuffle_spike_trains(simulated_spike_train(:,active_cells),num_shuffles,shuffle_type);

% Indetifying significantly modulated cells:
if measures_to_estimate(1) || measures_to_estimate(2) % based on the SI in active cells for naive versus shuffle
    % naive SI:
    [SI_naive_bit_spike,SI_naive_bit_sec]=compute_SI(average_firing_rates(active_cells),rate_maps(active_cells,:),normalized_states_distribution);
    
    % shuffle SI:
    SI_shuffle_bit_spike=nan(length(active_cells),num_shuffles);
    display_progress_bar('Computing shuffle information for the significance test: ',false)
    for n=1:num_shuffles
        display_progress_bar(100*(n/num_shuffles),false)
        [temp_shuffled_rate_maps,~,~]=compute_rate_maps(shuffled_spike_trains(:,:,n),stimulus_trace,dt);
        [SI_shuffle_bit_spike(:,n),~]=compute_SI(average_firing_rates(active_cells),temp_shuffled_rate_maps,normalized_states_distribution);
    end
    display_progress_bar(' done',false)
    display_progress_bar('',true)
    
    % Finding significant cells:
    information_significance_active_cells=1-sum(repmat(SI_naive_bit_spike,1,num_shuffles)>SI_shuffle_bit_spike,2)/num_shuffles;
    p_value_significant_active_cells=information_significance_active_cells(information_significance_active_cells<significance_threshold)';
    significant_active_cells_indexes=active_cells(information_significance_active_cells<significance_threshold);
    SI_naive_bit_spike_significant_cells=SI_naive_bit_spike(information_significance_active_cells<significance_threshold);
    SI_naive_bit_sec_significant_cells=SI_naive_bit_sec(information_significance_active_cells<significance_threshold);
    mean_SI_naive_bit_spike_significant_cells=mean(SI_naive_bit_spike_significant_cells,'omitnan');
    mean_SI_naive_bit_sec_significant_cells=mean(SI_naive_bit_sec_significant_cells,'omitnan');
    
elseif measures_to_estimate(3) % based on the MI in active cells for naive versus shuffle
    % naive MI:
    MI_naive=compute_MI(simulated_spike_train(:,active_cells),stimulus_trace);
    
    % shuffle MI:
    MI_shuffle=nan(length(active_cells),num_shuffles);
    display_progress_bar('Computing shuffle information: ',false)
    for n=1:num_shuffles
        display_progress_bar(100*(n/num_shuffles),false)
        MI_shuffle(:,n)=compute_MI(squeeze(shuffled_spike_trains(:,:,n)),stimulus_trace);
    end
    display_progress_bar(' done',false)
    display_progress_bar('',true)
    
    % Finding significant cells:
    information_significance_active_cells=1-sum(repmat(MI_naive,1,num_shuffles)>MI_shuffle,2)/num_shuffles;
    p_value_significant_active_cells=information_significance_active_cells(information_significance_active_cells<significance_threshold)';
    MI_naive_significant_cells=MI_naive(information_significance_active_cells<significance_threshold);
    mean_MI_naive_significant_cells=mean(MI_naive_significant_cells,'omitnan');
    significant_active_cells_indexes=active_cells(information_significance_active_cells<significance_threshold);
end

% Focusing only on significant cells:
disp(['Found ' num2str(length(active_cells)) '/' num2str(length(active_bins)) ' sufficiently active cells, out of which ' num2str(length(significant_active_cells_indexes)) ' cells (' num2str(round(100*length(significant_active_cells_indexes)/length(active_cells))) '%) are significantly modulated by the encoded variable.'])
average_rates_significant_cells=average_firing_rates(significant_active_cells_indexes);
active_bins_significant_cells=active_bins(significant_active_cells_indexes);
fraction_significant_and_active_cells=length(significant_active_cells_indexes)/length(active_bins);
fraction_of_significant_from_active_cells=length(significant_active_cells_indexes)/length(active_cells);
spike_trian_significant_cells=simulated_spike_train(:,significant_active_cells_indexes);
SI_true_bit_spike_significant_cells=true_information.SI_bit_spike(significant_active_cells_indexes);
SI_true_bit_sec_significant_cells=true_information.SI_bit_sec(significant_active_cells_indexes);
MI_true_significant_cells=true_information.MI(significant_active_cells_indexes);
mean_SI_true_bit_spike_significant_cells=mean(SI_true_bit_spike_significant_cells,'omitnan');
mean_SI_true_bit_sec_significant_cells=mean(SI_true_bit_sec_significant_cells,'omitnan');
mean_MI_true_significant_cells=mean(MI_true_significant_cells,'omitnan');

if (measures_to_estimate(1) || measures_to_estimate(2)) && measures_to_estimate(3) % compute naive MI even when not used for identifying significant cells
    MI_naive=compute_MI(simulated_spike_train(:,active_cells),stimulus_trace);
    MI_naive_significant_cells=MI_naive(information_significance_active_cells<significance_threshold);
    mean_MI_naive_significant_cells=mean(MI_naive_significant_cells,'omitnan');
end


%% Step 3: Calculating the naive and shuffle information as a function of sample size:

subsampling_repetitions=100; % number of repetitions in the subsampling of the data
T=size(simulated_spike_train,1); % total number of samples
subsample_size=(0.05:0.05:1)*T; % subsamples with size of different fractions of the data

disp('Computing information as a function of subsample size:')
if measures_to_estimate(1) || measures_to_estimate(2)
    if measures_to_estimate(3) % Compute SI and MI versus sample size
        [SI_naive_bit_spike_versus_sample_size,SI_shuffle_bit_spike_versus_sample_size,SI_naive_bit_sec_versus_sample_size,SI_shuffle_bit_sec_versus_sample_size,MI_naive_versus_sample_size,MI_shuffle_versus_sample_size]=...
            compute_information_versus_sample_size(spike_trian_significant_cells,stimulus_trace,subsample_size,dt,subsampling_repetitions,measures_to_estimate);
        SI_naive_bit_spike=SI_naive_bit_spike_versus_sample_size(:,end);
        SI_naive_bit_sec=SI_naive_bit_sec_versus_sample_size(:,end);
        MI_naive=MI_naive_versus_sample_size(:,end);
        
    else % Compute only SI versus sample size
        [SI_naive_bit_spike_versus_sample_size,SI_shuffle_bit_spike_versus_sample_size,SI_naive_bit_sec_versus_sample_size,SI_shuffle_bit_sec_versus_sample_size]=...
            compute_information_versus_sample_size(spike_trian_significant_cells,stimulus_trace,subsample_size,dt,subsampling_repetitions,measures_to_estimate);
    end
elseif measures_to_estimate(3) % Compute only MI versus sample size
    [~,~,~,~,MI_naive_versus_sample_size,MI_shuffle_versus_sample_size]=...
        compute_information_versus_sample_size(spike_trian_significant_cells,stimulus_trace,subsample_size,dt,subsampling_repetitions,measures_to_estimate);
end

%% Step 4: Correcting the bias  using the scaled shuffle reduction (SSR) and bounded extrapolation (BAE) methods:

% Thresholds for producing warnings:
SSR_stability_threshold=0.95;
BAE_fit_R_2_threshold=0.95;
mean_SSR_stability_threshold=0.98;
mean_BAE_fit_R_2_threshold=0.98;

% Correcting the bias for SI in bit/spike:
if measures_to_estimate(1)
    units='bit/spike';
    
    % SSR method:
    [SI_SSR_bit_spike,mean_SI_SSR_bit_spike,SI_SSR_stability_bit_spike,mean_SI_SSR_stability_bit_spike]...
        =perform_SSR_simulated_data(SI_naive_bit_spike_versus_sample_size,SI_shuffle_bit_spike_versus_sample_size,SI_true_bit_spike_significant_cells,subsample_size,units,SSR_stability_threshold,mean_SSR_stability_threshold,figures_directory);
    
    % BAE method:
    [SI_BAE_bit_spike,mean_SI_BAE_bit_spike,SI_BAE_fit_R_2_bit_spike,mean_SI_BAE_fit_R_2_bit_spike]...
        =perform_BAE_simulated_data(SI_naive_bit_spike_versus_sample_size,SI_true_bit_spike_significant_cells,subsample_size,units,BAE_fit_R_2_threshold,mean_BAE_fit_R_2_threshold,figures_directory);
    
    SI_disagreement_bit_spike=SI_BAE_bit_spike-SI_SSR_bit_spike;
    mean_SI_disagreement_bit_spike=mean_SI_BAE_bit_spike-mean_SI_SSR_bit_spike;
end

% Correcting the bias for SI in bit/sec:
if measures_to_estimate(2)
    units='bit/sec';
    
    % SSR method:
    [SI_SSR_bit_sec,mean_SI_SSR_bit_sec,SI_SSR_stability_bit_sec,mean_SI_SSR_stability_bit_sec]...
        =perform_SSR_simulated_data(SI_naive_bit_sec_versus_sample_size,SI_shuffle_bit_sec_versus_sample_size,SI_true_bit_sec_significant_cells,subsample_size,units,SSR_stability_threshold,mean_SSR_stability_threshold,figures_directory);
    
    % BAE method:
    [SI_BAE_bit_sec, mean_SI_BAE_bit_sec,SI_BAE_fit_R_2_bit_sec,mean_SI_BAE_fit_R_2_bit_sec]...
        =perform_BAE_simulated_data(SI_naive_bit_sec_versus_sample_size,SI_true_bit_sec_significant_cells,subsample_size,units,BAE_fit_R_2_threshold,mean_BAE_fit_R_2_threshold,figures_directory);
    
    SI_disagreement_bit_sec=SI_BAE_bit_sec-SI_SSR_bit_sec;
    mean_SI_disagreement_bit_sec=mean_SI_BAE_bit_sec-mean_SI_SSR_bit_sec;
end

% Correcting the bias for MI:
if measures_to_estimate(3)
    units='bit';
    
    % SSR method:
    [MI_SSR,mean_MI_SSR,MI_SSR_stability,mean_MI_SSR_stability]...
        =perform_SSR_simulated_data(MI_naive_versus_sample_size,MI_shuffle_versus_sample_size,MI_true_significant_cells,subsample_size,units,SSR_stability_threshold,mean_SSR_stability_threshold,figures_directory);
    
    % BAE method:
    [MI_BAE, mean_MI_BAE,MI_BAE_fit_R_2,mean_MI_BAE_fit_R_2]...
        =perform_BAE_simulated_data(MI_naive_versus_sample_size,MI_true_significant_cells,subsample_size,units,BAE_fit_R_2_threshold,mean_BAE_fit_R_2_threshold,figures_directory);
    
    MI_disagreement=MI_BAE-MI_SSR;
    mean_MI_disagreement=mean_MI_BAE-mean_MI_SSR;
end

%% Saving the final results in a single data structure:

% General parameters:
simulated_bias_correction_results=struct;
simulated_bias_correction_results.parameters=struct;
simulated_bias_correction_results.parameters.dt=dt;
simulated_bias_correction_results.parameters.num_shuffles=num_shuffles;
simulated_bias_correction_results.parameters.shuffle_type=shuffle_type;
simulated_bias_correction_results.parameters.active_bins_threshold=active_bins_threshold;
simulated_bias_correction_results.parameters.firing_rate_threshold=firing_rate_threshold;
simulated_bias_correction_results.parameters.significance_threshold=significance_threshold;
simulated_bias_correction_results.parameters.subsampling_repetitions=subsampling_repetitions;

% Data statistics:
simulated_bias_correction_results.statistics=struct;
simulated_bias_correction_results.statistics.fraction_significant_and_active_cells=fraction_significant_and_active_cells;
simulated_bias_correction_results.statistics.fraction_of_significant_from_active_cells=fraction_of_significant_from_active_cells;
simulated_bias_correction_results.statistics.average_rates_significant_cells=average_rates_significant_cells;
simulated_bias_correction_results.statistics.active_bins_significant_cells=active_bins_significant_cells;
simulated_bias_correction_results.statistics.significant_active_cells_indexes=significant_active_cells_indexes;
simulated_bias_correction_results.statistics.p_value_significant_active_cells=p_value_significant_active_cells;

% Estimated information:
simulated_bias_correction_results.information=struct;
if measures_to_estimate(1) % SI in bit/spike
    % For individual cells:
    simulated_bias_correction_results.information.SI_naive_bit_spike=SI_naive_bit_spike_significant_cells;
    simulated_bias_correction_results.information.SI_true_bit_spike=SI_true_bit_spike_significant_cells';
    simulated_bias_correction_results.information.SI_SSR_bit_spike=SI_SSR_bit_spike;
    simulated_bias_correction_results.information.SI_BAE_bit_spike=SI_BAE_bit_spike;
    simulated_bias_correction_results.information.SI_disagreement_bit_spike=SI_disagreement_bit_spike;
    simulated_bias_correction_results.information.SI_SSR_stability_bit_spike=SI_SSR_stability_bit_spike;
    simulated_bias_correction_results.information.SI_BAE_fit_R_2_bit_spike=SI_BAE_fit_R_2_bit_spike;
    
    % Average across the population:
    simulated_bias_correction_results.information.mean_SI_naive_bit_spike=mean_SI_naive_bit_spike_significant_cells;
    simulated_bias_correction_results.information.mean_SI_true_bit_spike=mean_SI_true_bit_spike_significant_cells;
    simulated_bias_correction_results.information.mean_SI_SSR_bit_spike=mean_SI_SSR_bit_spike;
    simulated_bias_correction_results.information.mean_SI_BAE_bit_spike=mean_SI_BAE_bit_spike;
    simulated_bias_correction_results.information.mean_SI_disagreement_bit_spike=mean_SI_disagreement_bit_spike;
    simulated_bias_correction_results.information.mean_SI_SSR_stability_bit_spike=mean_SI_SSR_stability_bit_spike;
    simulated_bias_correction_results.information.mean_SI_BAE_fit_R_2_bit_spike=mean_SI_BAE_fit_R_2_bit_spike;
end

if measures_to_estimate(2) % SI in bit/sec
    % For individual cells:
    simulated_bias_correction_results.information.SI_naive_bit_sec=SI_naive_bit_sec_significant_cells;
    simulated_bias_correction_results.information.SI_true_bit_sec=SI_true_bit_sec_significant_cells';
    simulated_bias_correction_results.information.SI_SSR_bit_sec=SI_SSR_bit_sec;
    simulated_bias_correction_results.information.SI_BAE_bit_sec=SI_BAE_bit_sec;
    simulated_bias_correction_results.information.SI_disagreement_bit_sec=SI_disagreement_bit_sec;
    simulated_bias_correction_results.information.SI_SSR_stability_bit_sec=SI_SSR_stability_bit_sec;
    simulated_bias_correction_results.information.SI_BAE_fit_R_2_bit_sec=SI_BAE_fit_R_2_bit_sec;
    
    % Average across the population:
    simulated_bias_correction_results.information.mean_SI_naive_bit_sec=mean_SI_naive_bit_sec_significant_cells;
    simulated_bias_correction_results.information.mean_SI_true_bit_sec=mean_SI_true_bit_sec_significant_cells;
    simulated_bias_correction_results.information.mean_SI_SSR_bit_sec=mean_SI_SSR_bit_sec;
    simulated_bias_correction_results.information.mean_SI_BAE_bit_sec=mean_SI_BAE_bit_sec;
    simulated_bias_correction_results.information.mean_SI_disagreement_bit_sec=mean_SI_disagreement_bit_sec;
    simulated_bias_correction_results.information.mean_SI_SSR_stability_bit_sec=mean_SI_SSR_stability_bit_sec;
    simulated_bias_correction_results.information.mean_SI_BAE_fit_R_2_bit_sec=mean_SI_BAE_fit_R_2_bit_sec;
end

if measures_to_estimate(3) % MI
    % For individual cells:
    simulated_bias_correction_results.information.MI_naive=MI_naive_significant_cells;
    simulated_bias_correction_results.information.MI_true=MI_true_significant_cells';
    simulated_bias_correction_results.information.MI_SSR=MI_SSR;
    simulated_bias_correction_results.information.MI_BAE=MI_BAE;
    simulated_bias_correction_results.information.MI_disagreement=MI_disagreement;
    simulated_bias_correction_results.information.MI_SSR_stability=MI_SSR_stability;
    simulated_bias_correction_results.information.MI_BAE_fit_R_2=MI_BAE_fit_R_2;
    
    % Average across the population:
    simulated_bias_correction_results.information.mean_MI_naive=mean_MI_naive_significant_cells;
    simulated_bias_correction_results.information.mean_MI_true=mean_MI_true_significant_cells;
    simulated_bias_correction_results.information.mean_MI_SSR=mean_MI_SSR;
    simulated_bias_correction_results.information.mean_MI_BAE=mean_MI_BAE;
    simulated_bias_correction_results.information.mean_MI_disagreement=mean_MI_disagreement;
    simulated_bias_correction_results.information.mean_MI_SSR_stability=mean_MI_SSR_stability;
    simulated_bias_correction_results.information.mean_MI_BAE_fit_R_2=mean_MI_BAE_fit_R_2;
end

save([data_path '\Sample data - results\simulated_bias_correction_results.mat'],'simulated_bias_correction_results');
disp('Finished analyzing data set')