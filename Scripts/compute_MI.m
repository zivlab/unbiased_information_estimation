function [MI]=compute_MI(spike_train,stimulus_trace)
% This function performs the naive calculation of the mutual information (MI) for a population of cells.

% Inputs:
% 1. spike_train - Matrix of size TxN, where T is the number of time bins and N is the number of neurons.
% Each element is the spike count of a given neuron in a given time bin.
% 2. stimulus_trace - Vector of size T, where each element is the state in a
% given time bin.

% Outputs:
% 1. MI - Vector of size N with the naive information of each cell

if isa(spike_train,'single') || isa(stimulus_trace,'single')
    epsilon=10^-30;
else
    epsilon=realmin;
end
states=unique(stimulus_trace);
num_states=length(states);
num_cells=size(spike_train,2);

% estimating the prior distribution of the encoded variable:
p_s=zeros(1,num_states);
for n=1:num_states
    p_s(n)=sum(stimulus_trace==states(n))/length(stimulus_trace);
end

% check if spike train consists of integer spike counts or a contiuous range of positive values:
is_integer_spike_count=sum(mod(spike_train(:),1))==0;
if is_integer_spike_count
    % binning the responses using integer values:
    max_rate=max(spike_train(:));
    response_bins=0:max_rate;
    
    % estimating the marginal probability of responses:
    p_r=zeros(num_cells,length(response_bins));
    for r=1:length(response_bins)
        p_r(:,r)=sum(spike_train==response_bins(r))/size(spike_train,1);
    end
    response_entropy=-sum(p_r.*log2(p_r+epsilon),2);
    
    % estimating the conditional probability of responses given a value of the encoded variable:
    p_r_given_s=zeros(num_cells,length(response_bins),num_states);
    for s=1:num_states
        for r=1:length(response_bins)
            if length(find(stimulus_trace==states(s)))>1
                p_r_given_s(:,r,s)=sum(spike_train(stimulus_trace==states(s),:)==response_bins(r))/length(find(stimulus_trace==states(s)));
            else
                p_r_given_s(:,r,s)=spike_train(stimulus_trace==states(s),:)==response_bins(r);
            end
        end
    end
else % binning the responses using non-integer values:
    if max(spike_train(:))<20 && max(spike_train(:))>5 % responses can be approximated by integer values
        spike_train=round(spike_train); % discretizing the spike trains using the closest integer values
        max_rate=max(spike_train(:));
        response_bins=0:max_rate;
    else % responses must binned using continuous ranges of values
        non_zero_responses=spike_train(spike_train>0);
        num_response_bins=11; % 1 bin for zero and 10 bins for the rest
        temp_non_zero_response_bins=logspace(log10(min(non_zero_responses)),log10(max(non_zero_responses)),num_response_bins);
        non_zero_response_bins=zeros(1,num_response_bins-1);
        for k=1:num_response_bins-1
            non_zero_response_bins(k)=mean(temp_non_zero_response_bins(k:k+1));
        end
        % discretizing the spike trains using the closest non-integer values:
        non_zero_spike_indexes=find(spike_train>0);
        for n=1:length(non_zero_spike_indexes)
           this_response_value=spike_train(non_zero_spike_indexes(n));
           [~,this_closest_response_bin]=min(abs(this_response_value-non_zero_response_bins));
           spike_train(non_zero_spike_indexes(n))=non_zero_response_bins(this_closest_response_bin);
        end
        response_bins=[0,non_zero_response_bins];
    end
    % estimating the marginal probability of responses:
    p_r=zeros(num_cells,length(response_bins));
    for r=1:length(response_bins)
        p_r(:,r)=sum(spike_train==response_bins(r))/size(spike_train,1);
    end
    response_entropy=-sum(p_r.*log2(p_r+epsilon),2);
    
    % estimating the conditional probability of responses given a value of the encoded variable:
    p_r_given_s=zeros(num_cells,length(response_bins),num_states);
    for s=1:num_states
        for r=1:length(response_bins)
            if length(find(stimulus_trace==states(s)))>1
                p_r_given_s(:,r,s)=sum(spike_train(stimulus_trace==states(s),:)==response_bins(r))/length(find(stimulus_trace==states(s)));
            else
                p_r_given_s(:,r,s)=spike_train(stimulus_trace==states(s),:)==response_bins(r);
            end
        end
    end
end

% computing the MI:
if num_cells>1
conditional_entropy=-sum(p_s.*squeeze(sum(p_r_given_s.*log2(p_r_given_s+epsilon),2)),2);
else
    conditional_entropy=-sum(p_s.*squeeze(sum(squeeze(p_r_given_s).*log2(squeeze(p_r_given_s)+epsilon),1)),2);
end
MI=response_entropy-conditional_entropy;
MI(isnan(MI))=0;

