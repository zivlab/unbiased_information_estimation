function shuffled_spike_trains=shuffle_spike_trains(spike_train,num_shuffles,shuffle_type)
% This function performs either a cyclic permutation or random shuffling 
% to obtain shuffled spike trains

% Inputs:
% 1. spike_train - Matrix of size TxN, where T is the number of time bins and N is the number of neurons. 
% Each element is the spike count of a given neuron in a given time bin.
% 2. num_shuffles - Number oif shuffling repetitions
% 3. shuffle_type - Either cyclic or random permutations

% Outputs:
% 1. shuffled_spike_trains - Matrix of size TxNxK, where T is the number of time bins, N is the number of neurons,
% and K is the number of shuffles. Each element is the spike count of a given neuron in a given time bin.

T=size(spike_train,1);
N=size(spike_train,2);
shuffled_spike_trains=zeros(T,N,num_shuffles);
if strcmp(shuffle_type,'cyclic')
    for n=1:num_shuffles
        shift_index=randi(T);
        shuffled_spike_trains(:,:,n)=spike_train([shift_index:end 1:shift_index-1],:);
    end
elseif strcmp(shuffle_type,'random')
    for n=1:num_shuffles
        [~,random_indexes]=sort(rand(1,T));
        shuffled_spike_trains(:,:,n)=spike_train(random_indexes,:);
    end
else
    warndlg('Please choose a valid shuffling type')
end

end
