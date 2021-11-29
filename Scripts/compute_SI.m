function [SI_bit_spike,SI_bit_sec]=compute_SI(average_firing_rates,tuning_curves,stimulus_distribution)
% This function performs the naive calculation of the Skaggs information index (SI) for a population of cells.

% Inputs: 
% 1. average_firing_rates - Vector of size N with the average firing rates of each cell
% 2. tuning_curves - Matrix of size NxS with the firing rate map of each
% 3. stimulus_distribution - Vector of size S with the probabilities of the different stimuli (stimuli)

% Outputs:
% 1. SI_bit_spike - Vector of size N with the naive Skaggs information in bit/spike of each cell
% 2. SI_bit_sec - Vector of size N with the naive Skaggs information in bit/sec of each cell

if isa(average_firing_rates,'single') || isa(tuning_curves,'single') || isa(stimulus_distribution,'single')
    epsilon=10^-30;
else
    epsilon=realmin;
end
normalized_tuning_curves=(tuning_curves)./(average_firing_rates+epsilon)';
SI_bit_spike=sum(stimulus_distribution.*normalized_tuning_curves.*log2(normalized_tuning_curves+epsilon),2);
SI_bit_spike(average_firing_rates==0)=nan;
SI_bit_sec=SI_bit_spike.*average_firing_rates';
SI_bit_sec(isnan(SI_bit_spike))=0;

end
