function [varargout]=compute_information_versus_sample_size(spike_train,stimulus_trace,subsample_size,dt,subsampling_repetitions,measures_to_estimate)
% This function estimates the population average and cell information
% as a function sample size using multiple subsamples of the data

% Inputs:
% 1. spike_train - Matrix of size TxN, where T is the number of time bins
% and N is the number of neurons. Each element is the spike count of a
% given neuron in a given time bin.
% 2. stimulus_trace - Vector of size T, where each element is the state in a
% given time bin.
% 3. subsample_size - Vector of T different sample sizes
% 4. dt - Temporal bin size (in seconds)
% 5. subsampling_repetitions - Number of subsampling repetitions for each
% sample duration.
% 6. measures_to_estimate - Vector of size 1x3 with values of 1 for computing the SI in bit/spike and
% bit/sec and the mutual information.

% Outputs:
% 1. varargout
% 1.1. SI_bit_spike_versus_sample_size - Matrix of size TxN with the estimated
% SI (bit/spike) of each of N neurons as a function of T different sample
% sizes.
% 1.2. SI_shuffle_bit_spike_versus_sample_size - Matrix of size TxN with the estimated
% shuffled SI (bit/spike) of each of N neurons as a function of T different
% sample sizes.
% 1.3. SI_bit_sec_versus_sample_size - Matrix of size TxN with the estimated
% SI (bit/sec) of each of N neurons as a function of T different sample
% sizes.
% 1.4. SI_shuffle_bit_sec_versus_sample_size - Matrix of size TxN with the estimated
% shuffled SI (bit/sec) of each of N neurons as a function of T different
% sample sizes.
% 1.5. MI_versus_sample_size - Matrix of size TxN with the estimated
% mutual information of each of N neurons as a function of T different sample
% sizes.
% 1.6. MI_shuffle_versus_sample_size - Matrix of size TxN with the estimated
% shuffled mutual information of each of N neurons as a function of T different
% sample sizes.

T=size(spike_train,1);
N=size(spike_train,2);

if measures_to_estimate(1) || measures_to_estimate(2)
    SI_naive_bit_spike_versus_sample_size=nan(N,length(subsample_size));
    SI_shuffle_bit_spike_versus_sample_size=nan(N,length(subsample_size));
    SI_naive_bit_sec_versus_sample_size=nan(N,length(subsample_size));
    SI_shuffle_bit_sec_versus_sample_size=nan(N,length(subsample_size));
end
if measures_to_estimate(3)
    MI_naive_versus_sample_size=nan(N,length(subsample_size));
    MI_shuffle_versus_sample_size=nan(N,length(subsample_size));
end

% calculating information across different subsample durations:
for n=1:length(subsample_size)
    this_num_time_bins=floor(subsample_size(n));
    if measures_to_estimate(1) ||  measures_to_estimate(2)
        this_subsample_size_SI_bit_spike=nan(N,ceil(subsampling_repetitions*T/subsample_size(n)));
        this_subsample_size_SI_shuffle_bit_spike=nan(N,ceil(subsampling_repetitions*T/subsample_size(n)));
        this_subsample_size_SI_bit_sec=nan(N,ceil(subsampling_repetitions*T/subsample_size(n)));
        this_subsample_size_SI_shuffle_bit_sec=nan(N,ceil(subsampling_repetitions*T/subsample_size(n)));
    end
    if measures_to_estimate(3)
        this_subsample_size_MI=nan(N,ceil(subsampling_repetitions*T/subsample_size(n)));
        this_subsample_size_MI_shuffle=nan(N,ceil(subsampling_repetitions*T/subsample_size(n)));
    end
    
    % Performing multiple repetitions of random subsampling for each duration:
    display_progress_bar('Terminating previous progress bars',true)
    display_progress_bar(['Computing information for a subsample size of ' num2str(n) '/' num2str(length(subsample_size)) ' of the full data size: '],false);
    for k=1:ceil(subsampling_repetitions*T/subsample_size(n))
        display_progress_bar(100*(k/ceil(subsampling_repetitions*T/subsample_size(n)) ),false)
        
        % Applying random shuffles to the data:
        [~,this_subsample_indexes]=sort(rand(T,1));
        this_subsample_indexes=this_subsample_indexes(1:this_num_time_bins);
        shuffled_spike_trains=shuffle_spike_trains(spike_train(this_subsample_indexes,:),1,'random');
        
        if measures_to_estimate(1) || measures_to_estimate(2)
            % computing SI for the subsampled rate maps:
            [temp_tuning_curves,temp_normalized_states_distribution]=compute_tuning_curves(spike_train(this_subsample_indexes,:),stimulus_trace(this_subsample_indexes),dt);
            temp_firing_rates=mean(spike_train(this_subsample_indexes,:))/dt;
            [temp_SI_bit_spike,temp_SI_bit_sec]=compute_SI(temp_firing_rates,temp_tuning_curves,temp_normalized_states_distribution);
            this_subsample_size_SI_bit_spike(:,k)=temp_SI_bit_spike;
            this_subsample_size_SI_bit_sec(:,k)=temp_SI_bit_sec;
            
            % computing SI for the shuffled subsample:
            [temp_shuffled_tuning_curves,~]=compute_tuning_curves(shuffled_spike_trains,stimulus_trace(this_subsample_indexes),dt);
            temp_shuffle_firing_rates=mean(shuffled_spike_trains)/dt;
            [temp_shuffle_information_bit_spike,temp_shuffle_information_bit_sec]=compute_SI(temp_shuffle_firing_rates,temp_shuffled_tuning_curves,temp_normalized_states_distribution);
            this_subsample_size_SI_shuffle_bit_spike(:,k)=temp_shuffle_information_bit_spike;
            this_subsample_size_SI_shuffle_bit_sec(:,k)=temp_shuffle_information_bit_sec;
        end
        if measures_to_estimate(3)
            % computing mutual information for the subsampled spike trains:
            temp_MI=compute_MI(spike_train(this_subsample_indexes,:),stimulus_trace(this_subsample_indexes));
            this_subsample_size_MI(:,k)=temp_MI;
            
            % computing mutual information for the shuffled subsampled spike train:
            temp_MI_shuffle=compute_MI(shuffled_spike_trains,stimulus_trace(this_subsample_indexes));
            this_subsample_size_MI_shuffle(:,k)=temp_MI_shuffle;
        end
    end
    
    % averaging the information across all subsample repetitions:
    if measures_to_estimate(1) || measures_to_estimate(2)
        SI_naive_bit_spike_versus_sample_size(:,n)=mean(this_subsample_size_SI_bit_spike,2,'omitnan');
        SI_shuffle_bit_spike_versus_sample_size(:,n)=mean(this_subsample_size_SI_shuffle_bit_spike,2,'omitnan');
        SI_naive_bit_sec_versus_sample_size(:,n)=mean(this_subsample_size_SI_bit_sec,2,'omitnan');
        SI_shuffle_bit_sec_versus_sample_size(:,n)=mean(this_subsample_size_SI_shuffle_bit_sec,2,'omitnan');
    end
    if measures_to_estimate(3)
        MI_naive_versus_sample_size(:,n)=mean(this_subsample_size_MI,2,'omitnan');
        MI_shuffle_versus_sample_size(:,n)=mean(this_subsample_size_MI_shuffle,2,'omitnan');
    end
    display_progress_bar(' done',false)
    display_progress_bar('',true)
end

if measures_to_estimate(1) || measures_to_estimate(2)
    varargout{1}=SI_naive_bit_spike_versus_sample_size;
    varargout{2}=SI_shuffle_bit_spike_versus_sample_size;
    varargout{3}=SI_naive_bit_sec_versus_sample_size;
    varargout{4}=SI_shuffle_bit_sec_versus_sample_size;
end
if measures_to_estimate(3)
    varargout{5}=MI_naive_versus_sample_size;
    varargout{6}=MI_shuffle_versus_sample_size;
end