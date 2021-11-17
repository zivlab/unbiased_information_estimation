function unbiased_information_estimation_results=estimate_unbiased_information(spike_train,stimulus_trace,settings)
% This function performs all the steps for obtaining a bias-free estimation
% of information content.

% Inputs:
% 1. spike_train - Matrix of size TxN, where T is the number of time bins and N is the number of neurons.
% Each element is the spike count of a given neuron in a given time bin.
% 2. stimulus_trace - Vector of size T, where each element is the state in a
% given time bin.
% 3. settings - all the settings and parameters for performing the estimation:

% Outputs:
% 1. unbiased_information_estimation_results - all the results in a single data structure

measures_to_estimate=settings.measures_to_estimate;
dt=settings.dt;
active_bins_threshold=settings.active_bins_threshold;
firing_rate_threshold=settings.firing_rate_threshold;
active_bins=sum(spike_train>0);
average_firing_rates=mean(spike_train)/dt;
active_cells=find(active_bins>=active_bins_threshold & average_firing_rates>firing_rate_threshold);

% Computing the naive SI:
if measures_to_estimate(1) || measures_to_estimate(2)
    % Computing the tuning curves of the cells:
    [tuning_curves,normalized_states_distribution]=compute_tuning_curves(spike_train,stimulus_trace,dt);
    % naive SI:
    [SI_naive_bit_spike,SI_naive_bit_sec]=compute_SI(average_firing_rates(active_cells),tuning_curves(active_cells,:),normalized_states_distribution);
end

% Computing the naive MI:
if measures_to_estimate(3) % compute naive MI even when not used for identifying significant cells
    MI_naive=compute_MI(spike_train(:,active_cells),stimulus_trace);
end

% Identifying significantly modulated cells by comparing the naive information with shuffles:
estimate_only_significant_cells=settings.estimate_only_significant_cells;
if estimate_only_significant_cells
    shuffle_type=settings.shuffle_type;
    num_shuffles=settings.num_shuffles;
    significance_threshold=settings.significance_threshold;
    
    % obtaining shuffled spike trains:
    shuffled_spike_trains=shuffle_spike_trains(spike_train(:,active_cells),num_shuffles,shuffle_type);
    
    % Indetifying significantly modulated cells:
    if measures_to_estimate(1) || measures_to_estimate(2) % based on the SI in active cells for naive versus shuffle
        
        % shuffle SI:
        SI_shuffle_bit_spike=nan(length(active_cells),num_shuffles);
        display_progress_bar('Computing shuffle information for the significance test: ',false)
        for n=1:num_shuffles
            display_progress_bar(100*(n/num_shuffles),false)
            [temp_shuffled_tuning_curves,~]=compute_tuning_curves(shuffled_spike_trains(:,:,n),stimulus_trace,dt);
            [SI_shuffle_bit_spike(:,n),~]=compute_SI(average_firing_rates(active_cells),temp_shuffled_tuning_curves,normalized_states_distribution);
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
    
    if (measures_to_estimate(1) || measures_to_estimate(2)) && measures_to_estimate(3) % compute naive MI even when not used for identifying significant cells
        MI_naive_significant_cells=MI_naive(information_significance_active_cells<significance_threshold);
        mean_MI_naive_significant_cells=mean(MI_naive_significant_cells,'omitnan');
    end
else
    significant_active_cells_indexes=active_cells;
    if measures_to_estimate(1) || measures_to_estimate(2)
        SI_naive_bit_spike_significant_cells=SI_naive_bit_spike;
        SI_naive_bit_sec_significant_cells=SI_naive_bit_sec;
        mean_SI_naive_bit_spike_significant_cells=mean(SI_naive_bit_spike_significant_cells,'omitnan');
        mean_SI_naive_bit_sec_significant_cells=mean(SI_naive_bit_sec_significant_cells,'omitnan');
    elseif measures_to_estimate(3)
        MI_naive_significant_cells=MI_naive;
        mean_MI_naive_significant_cells=mean(MI_naive_significant_cells,'omitnan');
    end
    if (measures_to_estimate(1) || measures_to_estimate(2)) && measures_to_estimate(3) % compute naive MI even when not used for identifying significant cells
        MI_naive_significant_cells=MI_naive;
        mean_MI_naive_significant_cells=mean(MI_naive_significant_cells,'omitnan');
    end
end

% Focusing only on significant cells:
average_rates_significant_cells=average_firing_rates(significant_active_cells_indexes);
active_bins_significant_cells=active_bins(significant_active_cells_indexes);
fraction_significant_and_active_cells=length(significant_active_cells_indexes)/length(active_bins);
fraction_of_significant_from_active_cells=length(significant_active_cells_indexes)/length(active_cells);
spike_train_significant_cells=spike_train(:,significant_active_cells_indexes);
if estimate_only_significant_cells
    disp(['Found ' num2str(length(active_cells)) '/' num2str(length(active_bins)) ' sufficiently active cells, out of which ' num2str(length(significant_active_cells_indexes)) ' cells (' num2str(round(100*length(significant_active_cells_indexes)/length(active_cells))) '%) are significantly modulated by the encoded variable.'])
else
    disp(['Found ' num2str(length(active_cells)) '/' num2str(length(active_bins)) ' sufficiently active cells'])
end

% Calculating the naive and shuffle information as a function of sample size:
subsampling_repetitions=settings.subsampling_repetitions; % number of repetitions in the subsampling of the data
T=size(spike_train,1); % total number of samples
subsample_size=settings.subsample_fraction*T; % subsamples with size of different fractions of the data

disp('Computing information as a function of subsample size:')
if measures_to_estimate(1) || measures_to_estimate(2)
    if measures_to_estimate(3) % Compute SI and MI versus sample size
        [SI_naive_bit_spike_versus_sample_size,SI_shuffle_bit_spike_versus_sample_size,SI_naive_bit_sec_versus_sample_size,SI_shuffle_bit_sec_versus_sample_size,MI_naive_versus_sample_size,MI_shuffle_versus_sample_size]=...
            compute_information_versus_sample_size(spike_train_significant_cells,stimulus_trace,subsample_size,dt,subsampling_repetitions,measures_to_estimate);     
    else % Compute only SI versus sample size
        [SI_naive_bit_spike_versus_sample_size,SI_shuffle_bit_spike_versus_sample_size,SI_naive_bit_sec_versus_sample_size,SI_shuffle_bit_sec_versus_sample_size]=...
            compute_information_versus_sample_size(spike_train_significant_cells,stimulus_trace,subsample_size,dt,subsampling_repetitions,measures_to_estimate);
    end
elseif measures_to_estimate(3) % Compute only MI versus sample size
    [~,~,~,~,MI_naive_versus_sample_size,MI_shuffle_versus_sample_size]=...
        compute_information_versus_sample_size(spike_train_significant_cells,stimulus_trace,subsample_size,dt,subsampling_repetitions,measures_to_estimate);
end

% Correcting the bias using the SSR BAE methods:
plot_results=settings.plot_results;
figures_directory=settings.figures_directory;
if measures_to_estimate(1) % Correcting the bias for SI in bit/spike

    units='bit/spike';  
    
    % SSR method:
    [SI_SSR_bit_spike,mean_SI_SSR_bit_spike,SI_SSR_stability_bit_spike,mean_SI_SSR_stability_bit_spike]...
        =perform_SSR(SI_naive_bit_spike_versus_sample_size,SI_shuffle_bit_spike_versus_sample_size,subsample_size,units,figures_directory,plot_results);
    
    % BAE method:
    [SI_BAE_bit_spike,mean_SI_BAE_bit_spike,SI_BAE_fit_R_2_bit_spike,mean_SI_BAE_fit_R_2_bit_spike]...
        =perform_BAE(SI_naive_bit_spike_versus_sample_size,subsample_size,units,figures_directory,plot_results);
    
    SI_disagreement_bit_spike=SI_BAE_bit_spike-SI_SSR_bit_spike;
    mean_SI_disagreement_bit_spike=mean_SI_BAE_bit_spike-mean_SI_SSR_bit_spike;
end

% Correcting the bias for SI in bit/sec:
if measures_to_estimate(2)
    units='bit/sec';
    
    % SSR method:
    [SI_SSR_bit_sec,mean_SI_SSR_bit_sec,SI_SSR_stability_bit_sec,mean_SI_SSR_stability_bit_sec]...
        =perform_SSR(SI_naive_bit_sec_versus_sample_size,SI_shuffle_bit_sec_versus_sample_size,subsample_size,units,figures_directory,plot_results);
    
    % BAE method:
    [SI_BAE_bit_sec, mean_SI_BAE_bit_sec,SI_BAE_fit_R_2_bit_sec,mean_SI_BAE_fit_R_2_bit_sec]...
        =perform_BAE(SI_naive_bit_sec_versus_sample_size,subsample_size,units,figures_directory,plot_results);
    
    SI_disagreement_bit_sec=SI_BAE_bit_sec-SI_SSR_bit_sec;
    mean_SI_disagreement_bit_sec=mean_SI_BAE_bit_sec-mean_SI_SSR_bit_sec;
end

% Correcting the bias for MI:
if measures_to_estimate(3)
    units='bit';
    
    % SSR method:
    [MI_SSR,mean_MI_SSR,MI_SSR_stability,mean_MI_SSR_stability]...
        =perform_SSR(MI_naive_versus_sample_size,MI_shuffle_versus_sample_size,subsample_size,units,figures_directory,plot_results);
    
    % BAE method:
    [MI_BAE, mean_MI_BAE,MI_BAE_fit_R_2,mean_MI_BAE_fit_R_2]...
        =perform_BAE(MI_naive_versus_sample_size,subsample_size,units,figures_directory,plot_results);
    
    MI_disagreement=MI_BAE-MI_SSR;
    mean_MI_disagreement=mean_MI_BAE-mean_MI_SSR;
end

% plotting the results:
if plot_results
    if measures_to_estimate(1) % for SI in bit/spike
        
        % Cross validation:
        figure
        plot(SI_SSR_bit_spike,SI_BAE_bit_spike,'.','markersize',15,'color','b')
        hold on
        plot([0 1.1*max(SI_BAE_bit_spike)],[0 1.1*max(SI_BAE_bit_spike)],'--k','linewidth',2)
        xlim([0 1.1*max(SI_BAE_bit_spike)])
        ylim([0 1.1*max(SI_BAE_bit_spike)])
        axis square
        xlabel('SSR estimation (bit/spike)')
        ylabel('BAE estimation (bit/spike)')
        title('SI (bit/spike)')
        set(gca,'fontsize',16)
        box off
        savefig(fullfile(figures_directory,'BAE versus SSR - SI bit per spike.fig'))
        saveas(gcf,fullfile(figures_directory,'BAE versus SSR - SI bit per spike'),'png')
        
        % Estimation quality
        figure
        plot(SI_SSR_stability_bit_spike,SI_BAE_fit_R_2_bit_spike,'.','markersize',15,'color','b')
        xlim([0 1])
        ylim([0 1])
        axis square
        xlabel('SSR stability')
        ylabel('BAE fit R^2')
        title('SI (bit/spike)')
        set(gca,'fontsize',16)
        box off
        savefig(fullfile(figures_directory,'Estimation quality - SI bit per spike.fig'))
        saveas(gcf,fullfile(figures_directory,'Estimation quality - SI bit per spike'),'png')
    end
    
    if measures_to_estimate(2) % for SI in bit/sec
        
        % Cross validation:
        figure
        plot(SI_SSR_bit_sec,SI_BAE_bit_sec,'.','markersize',15,'color','b')
        hold on
        plot([0 1.1*max(SI_BAE_bit_sec)],[0 1.1*max(SI_BAE_bit_sec)],'--k','linewidth',2)
        xlim([0 1.1*max(SI_BAE_bit_sec)])
        ylim([0 1.1*max(SI_BAE_bit_sec)])
        axis square
        xlabel('SSR estimation (bit/sec)')
        ylabel('BAE estimation (bit/sec)')
        title('SI (bit/sec)')
        set(gca,'fontsize',16)
        box off
        savefig(fullfile(figures_directory,'BAE versus SSR - SI bit per sec.fig'))
        saveas(gcf,fullfile(figures_directory,'BAE versus SSR - SI bit per sec'),'png')
        
        % Estimation quality
        figure
        plot(SI_SSR_stability_bit_sec,SI_BAE_fit_R_2_bit_sec,'.','markersize',15,'color','b')
        xlim([0 1])
        ylim([0 1])
        axis square
        xlabel('SSR stability')
        ylabel('BAE fit R^2')
        title('SI (bit/sec)')
        set(gca,'fontsize',16)
        box off
        savefig(fullfile(figures_directory,'Estimation quality - SI bit per sec.fig'))
        saveas(gcf,fullfile(figures_directory,'Estimation quality - SI bit per sec'),'png')
    end
    
    if measures_to_estimate(3) % for MI
        
        % Cross validation:
        figure
        plot(MI_SSR,MI_BAE,'.','markersize',15,'color','b')
        hold on
        plot([0 1.1*max(MI_BAE)],[0 1.1*max(MI_BAE)],'--k','linewidth',2)
        xlim([0 1.1*max(MI_BAE)])
        ylim([0 1.1*max(MI_BAE)])
        axis square
        xlabel('SSR estimation (bit)')
        ylabel('BAE estimation (bit)')
        title('MI')
        set(gca,'fontsize',16)
        box off
        savefig(fullfile(figures_directory,'BAE versus SSR - MI.fig'))
        saveas(gcf,fullfile(figures_directory,'BAE versus SSR - MI'),'png')
        
        % Estimation quality
        figure
        plot(MI_SSR_stability,MI_BAE_fit_R_2,'.','markersize',15,'color','b')
        xlim([0 1])
        ylim([0 1])
        axis square
        xlabel('SSR stability')
        ylabel('BAE fit R^2')
        title('MI')
        set(gca,'fontsize',16)
        box off
        savefig(fullfile(figures_directory,'Estimation quality - MI.fig'))
        saveas(gcf,fullfile(figures_directory,'Estimation quality - MI'),'png')
    end
end

% Saving the final results in a single data structure:
% General parameters:
unbiased_information_estimation_results=struct;
unbiased_information_estimation_results.parameters=struct;
unbiased_information_estimation_results.parameters.dt=dt;
unbiased_information_estimation_results.parameters.num_shuffles=num_shuffles;
unbiased_information_estimation_results.parameters.shuffle_type=shuffle_type;
unbiased_information_estimation_results.parameters.active_bins_threshold=active_bins_threshold;
unbiased_information_estimation_results.parameters.firing_rate_threshold=firing_rate_threshold;
unbiased_information_estimation_results.parameters.significance_threshold=significance_threshold;
unbiased_information_estimation_results.parameters.subsampling_repetitions=subsampling_repetitions;

% Data statistics:
unbiased_information_estimation_results.statistics=struct;
unbiased_information_estimation_results.statistics.fraction_significant_and_active_cells=fraction_significant_and_active_cells;
unbiased_information_estimation_results.statistics.fraction_of_significant_from_active_cells=fraction_of_significant_from_active_cells;
unbiased_information_estimation_results.statistics.average_rates_significant_cells=average_rates_significant_cells;
unbiased_information_estimation_results.statistics.active_bins_significant_cells=active_bins_significant_cells;
unbiased_information_estimation_results.statistics.significant_active_cells_indexes=significant_active_cells_indexes;
unbiased_information_estimation_results.statistics.p_value_significant_active_cells=p_value_significant_active_cells;

% Estimated information:
unbiased_information_estimation_results.information=struct;
if measures_to_estimate(1) % SI in bit/spike
    % For individual cells:
    unbiased_information_estimation_results.information.SI_naive_bit_spike=SI_naive_bit_spike_significant_cells;
    unbiased_information_estimation_results.information.SI_SSR_bit_spike=SI_SSR_bit_spike;
    unbiased_information_estimation_results.information.SI_BAE_bit_spike=SI_BAE_bit_spike;
    unbiased_information_estimation_results.information.SI_disagreement_bit_spike=SI_disagreement_bit_spike;
    unbiased_information_estimation_results.information.SI_SSR_stability_bit_spike=SI_SSR_stability_bit_spike;
    unbiased_information_estimation_results.information.SI_BAE_fit_R_2_bit_spike=SI_BAE_fit_R_2_bit_spike;
    
    % Average across the population:
    unbiased_information_estimation_results.information.mean_SI_naive_bit_spike=mean_SI_naive_bit_spike_significant_cells;
    unbiased_information_estimation_results.information.mean_SI_SSR_bit_spike=mean_SI_SSR_bit_spike;
    unbiased_information_estimation_results.information.mean_SI_BAE_bit_spike=mean_SI_BAE_bit_spike;
    unbiased_information_estimation_results.information.mean_SI_disagreement_bit_spike=mean_SI_disagreement_bit_spike;
    unbiased_information_estimation_results.information.mean_SI_SSR_stability_bit_spike=mean_SI_SSR_stability_bit_spike;
    unbiased_information_estimation_results.information.mean_SI_BAE_fit_R_2_bit_spike=mean_SI_BAE_fit_R_2_bit_spike;
    
    % Average across spikes (weighted by the cells firing rates):
    unbiased_information_estimation_results.information.weighted_mean_SI_naive_bit_spike=sum(SI_naive_bit_spike_significant_cells.*average_rates_significant_cells')/sum(average_rates_significant_cells);
    unbiased_information_estimation_results.information.weighted_mean_SI_SSR_bit_spike=sum(SI_SSR_bit_spike.*average_rates_significant_cells')/sum(average_rates_significant_cells);
    unbiased_information_estimation_results.information.weighted_mean_SI_BAE_bit_spike=sum(SI_BAE_bit_spike.*average_rates_significant_cells')/sum(average_rates_significant_cells);
    unbiased_information_estimation_results.information.weighted_mean_SI_disagreement_bit_spike=sum(SI_disagreement_bit_spike.*average_rates_significant_cells')/sum(average_rates_significant_cells);
end

if measures_to_estimate(2) % SI in bit/sec
    % For individual cells:
    unbiased_information_estimation_results.information.SI_naive_bit_sec=SI_naive_bit_sec_significant_cells;
    unbiased_information_estimation_results.information.SI_SSR_bit_sec=SI_SSR_bit_sec;
    unbiased_information_estimation_results.information.SI_BAE_bit_sec=SI_BAE_bit_sec;
    unbiased_information_estimation_results.information.SI_disagreement_bit_sec=SI_disagreement_bit_sec;
    unbiased_information_estimation_results.information.SI_SSR_stability_bit_sec=SI_SSR_stability_bit_sec;
    unbiased_information_estimation_results.information.SI_BAE_fit_R_2_bit_sec=SI_BAE_fit_R_2_bit_sec;
    
    % Average across the population:
    unbiased_information_estimation_results.information.mean_SI_naive_bit_sec=mean_SI_naive_bit_sec_significant_cells;
    unbiased_information_estimation_results.information.mean_SI_SSR_bit_sec=mean_SI_SSR_bit_sec;
    unbiased_information_estimation_results.information.mean_SI_BAE_bit_sec=mean_SI_BAE_bit_sec;
    unbiased_information_estimation_results.information.mean_SI_disagreement_bit_sec=mean_SI_disagreement_bit_sec;
    unbiased_information_estimation_results.information.mean_SI_SSR_stability_bit_sec=mean_SI_SSR_stability_bit_sec;
    unbiased_information_estimation_results.information.mean_SI_BAE_fit_R_2_bit_sec=mean_SI_BAE_fit_R_2_bit_sec;
end

if measures_to_estimate(3) % MI
    % For individual cells:
    unbiased_information_estimation_results.information.MI_naive=MI_naive_significant_cells;
    unbiased_information_estimation_results.information.MI_SSR=MI_SSR;
    unbiased_information_estimation_results.information.MI_BAE_=MI_BAE;
    unbiased_information_estimation_results.information.MI_disagreement=MI_disagreement;
    unbiased_information_estimation_results.information.MI_SSR_stability=MI_SSR_stability;
    unbiased_information_estimation_results.information.MI_BAE_fit_R_2=MI_BAE_fit_R_2;
    
    % Average across the population:
    unbiased_information_estimation_results.information.mean_MI_naive=mean_MI_naive_significant_cells;
    unbiased_information_estimation_results.information.mean_MI_SSR=mean_MI_SSR;
    unbiased_information_estimation_results.information.mean_MI_BAE_=mean_MI_BAE;
    unbiased_information_estimation_results.information.mean_MI_disagreement=mean_MI_disagreement;
    unbiased_information_estimation_results.information.mean_MI_SSR_stability=mean_MI_SSR_stability;
    unbiased_information_estimation_results.information.mean_MI_BAE_fit_R_2=mean_MI_BAE_fit_R_2;
end

end

