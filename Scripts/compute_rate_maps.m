function [rate_maps,average_firing_rates,stimulus_distribution]=compute_rate_maps(spike_train,stimulus_trace,dt)

% Inputs: 
% 1. spike_train - Matrix of size TxN, where T is the number of time bins
% and N is the number of neurons. Each element is the spike count of a
% given neuron in a given time bin.
% 2. stimulus_trace - Vector of size T, where each element is the state in a
% 3. dt - Temporal bin size (in seconds) 

% Outputs:
% 1. average_firing_rates - Vector of size N with the average firing rate of each
% neuron 
% 2. rate_maps - Matrix of size NxS with the firing rate map of each
% 3. stimulus_distribution - Vector of size S with the probabilities of the different stimuli

stimulus_values=unique(stimulus_trace);
num_stimulus_values=length(stimulus_values);
num_cells=size(spike_train,2);

stimulus_distribution=hist(stimulus_trace,stimulus_values)/length(stimulus_trace);
average_firing_rates=mean(spike_train)/dt;
rate_maps=zeros(num_cells,num_stimulus_values);
for n=1:num_stimulus_values
    this_bin_indexes=find(stimulus_trace==stimulus_values(n));    
    if length(this_bin_indexes)>1
        rate_maps(:,n)=mean(spike_train(this_bin_indexes,:))/dt;
    else
        rate_maps(:,n)=spike_train(this_bin_indexes,:)/dt;
    end
end

end

