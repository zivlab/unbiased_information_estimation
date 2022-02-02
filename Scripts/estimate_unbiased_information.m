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

unbiased_information_estimation_results=struct;

measures_to_estimate=settings.measures_to_estimate;
dt=settings.dt;
active_bins_threshold=settings.active_bins_threshold;
firing_rate_threshold=settings.firing_rate_threshold;
active_bins=sum(spike_train>0);
average_firing_rates=mean(spike_train)/dt;
sufficiently_active_cells_indexes=find(active_bins>=active_bins_threshold & average_firing_rates>firing_rate_threshold);
fraction_sufficiently_active_cells=length(sufficiently_active_cells_indexes)/length(active_bins);

% Computing the naive SI:
if ~isempty(sufficiently_active_cells_indexes)
    if measures_to_estimate(1) || measures_to_estimate(2)
        % Computing the tuning curves of the cells:
        [tuning_curves,normalized_states_distribution]=compute_tuning_curves(spike_train,stimulus_trace,dt);
        % naive SI:
        [SI_naive_bit_spike,SI_naive_bit_sec]=compute_SI(average_firing_rates(sufficiently_active_cells_indexes),tuning_curves(sufficiently_active_cells_indexes,:),normalized_states_distribution);
    end
    
    % Computing the naive MI:
    if measures_to_estimate(3) % compute naive MI even when not used for identifying significantly tuned cells
        MI_naive=compute_MI(spike_train(:,sufficiently_active_cells_indexes),stimulus_trace);
    end
    
    % Identifying significantly modulated cells by comparing the naive information with shuffles:
    display_progress_bar('Terminating previous progress bars',true)
    estimate_only_significantly_tuned_cells=settings.estimate_only_significantly_tuned_cells;
    if estimate_only_significantly_tuned_cells
        shuffle_type=settings.shuffle_type;
        num_shuffles=settings.num_shuffles;
        tuning_significance_threshold=settings.tuning_significance_threshold;
        
        % obtaining shuffled spike trains:
        shuffled_spike_trains=shuffle_spike_trains(spike_train(:,sufficiently_active_cells_indexes),num_shuffles,shuffle_type);
        
        % Indetifying significantly modulated cells:
        if measures_to_estimate(1) || measures_to_estimate(2) % based on the SI in active cells for naive versus shuffle
            
            % shuffle SI:
            SI_shuffle_bit_spike=nan(length(sufficiently_active_cells_indexes),num_shuffles);
            display_progress_bar('Computing shuffle information for the tuning significance test: ',false)
            for n=1:num_shuffles
                display_progress_bar(100*(n/num_shuffles),false)
                [temp_shuffled_tuning_curves,~]=compute_tuning_curves(shuffled_spike_trains(:,:,n),stimulus_trace,dt);
                [SI_shuffle_bit_spike(:,n),~]=compute_SI(average_firing_rates(sufficiently_active_cells_indexes),temp_shuffled_tuning_curves,normalized_states_distribution);
            end
            display_progress_bar(' done',false)
            display_progress_bar('',true)
            
            % Finding significantly tuned cells:
            tuning_significance_active_cells=1-sum(repmat(SI_naive_bit_spike,1,num_shuffles)>SI_shuffle_bit_spike,2)/num_shuffles;
            p_value_significantly_tuned_and_active_cells=tuning_significance_active_cells(tuning_significance_active_cells<tuning_significance_threshold)';
            significantly_tuned_and_active_cells_indexes=sufficiently_active_cells_indexes(tuning_significance_active_cells<tuning_significance_threshold);
            SI_naive_bit_spike_significantly_tuned_cells=SI_naive_bit_spike(tuning_significance_active_cells<tuning_significance_threshold);
            SI_naive_bit_sec_significantly_tuned_cells=SI_naive_bit_sec(tuning_significance_active_cells<tuning_significance_threshold);
            if ~isempty(significantly_tuned_and_active_cells_indexes)
                average_SI_naive_bit_spike_significantly_tuned_cells=mean(SI_naive_bit_spike_significantly_tuned_cells,'omitnan');
                average_SI_naive_bit_sec_significantly_tuned_cells=mean(SI_naive_bit_sec_significantly_tuned_cells,'omitnan');
            end
        elseif measures_to_estimate(3) % based on the MI in active cells for naive versus shuffle
            
            % shuffle MI:
            MI_shuffle=nan(length(sufficiently_active_cells_indexes),num_shuffles);
            display_progress_bar('Computing shuffle information: ',false)
            for n=1:num_shuffles
                display_progress_bar(100*(n/num_shuffles),false)
                MI_shuffle(:,n)=compute_MI(squeeze(shuffled_spike_trains(:,:,n)),stimulus_trace);
            end
            display_progress_bar(' done',false)
            display_progress_bar('',true)
            
            % Finding significantly tuned cells:
            tuning_significance_active_cells=1-sum(repmat(MI_naive,1,num_shuffles)>MI_shuffle,2)/num_shuffles;
            p_value_significantly_tuned_and_active_cells=tuning_significance_active_cells(tuning_significance_active_cells<tuning_significance_threshold)';
            MI_naive_significantly_tuned_cells=MI_naive(tuning_significance_active_cells<tuning_significance_threshold);
            if ~isempty(MI_naive_significantly_tuned_cells)
                average_MI_naive_significantly_tuned_cells=mean(MI_naive_significantly_tuned_cells,'omitnan');
                significantly_tuned_and_active_cells_indexes=sufficiently_active_cells_indexes(tuning_significance_active_cells<tuning_significance_threshold);
            end
        end
        
        if (measures_to_estimate(1) || measures_to_estimate(2)) && measures_to_estimate(3) % compute naive MI even when not used for identifying significantly tuned cells
            MI_naive_significantly_tuned_cells=MI_naive(tuning_significance_active_cells<tuning_significance_threshold);
            if ~isempty(MI_naive_significantly_tuned_cells)
                average_MI_naive_significantly_tuned_cells=mean(MI_naive_significantly_tuned_cells,'omitnan');
            end
        end
    else
        significantly_tuned_and_active_cells_indexes=sufficiently_active_cells_indexes;
        if measures_to_estimate(1) || measures_to_estimate(2)
            SI_naive_bit_spike_significantly_tuned_cells=SI_naive_bit_spike;
            SI_naive_bit_sec_significantly_tuned_cells=SI_naive_bit_sec;
            if ~isempty(SI_naive_bit_spike_significantly_tuned_cells)
                average_SI_naive_bit_spike_significantly_tuned_cells=mean(SI_naive_bit_spike_significantly_tuned_cells,'omitnan');
                average_SI_naive_bit_sec_significantly_tuned_cells=mean(SI_naive_bit_sec_significantly_tuned_cells,'omitnan');
            end
        elseif measures_to_estimate(3)
            MI_naive_significantly_tuned_cells=MI_naive;
            if ~isempty(MI_naive_significantly_tuned_cells)
                average_MI_naive_significantly_tuned_cells=mean(MI_naive_significantly_tuned_cells,'omitnan');
            end
        end
        if (measures_to_estimate(1) || measures_to_estimate(2)) && measures_to_estimate(3) % compute naive MI even when not used for identifying significantly tuned cells
            MI_naive_significantly_tuned_cells=MI_naive;
            if ~isempty(MI_naive_significantly_tuned_cells)
                average_MI_naive_significantly_tuned_cells=mean(MI_naive_significantly_tuned_cells,'omitnan');
            end
        end
    end
    
    % Focusing only on significantly tuned cells:
    if ~isempty(significantly_tuned_and_active_cells_indexes)
        average_rates_significantly_tuned_cells=average_firing_rates(significantly_tuned_and_active_cells_indexes);
        active_bins_significantly_tuned_cells=active_bins(significantly_tuned_and_active_cells_indexes);
        fraction_significantly_tuned_and_active_cells=length(significantly_tuned_and_active_cells_indexes)/length(active_bins);
        fraction_of_significantly_tuned_from_active_cells=length(significantly_tuned_and_active_cells_indexes)/length(sufficiently_active_cells_indexes);
        spike_train_significantly_tuned_cells=spike_train(:,significantly_tuned_and_active_cells_indexes);
        if estimate_only_significantly_tuned_cells
            disp(['Found ' num2str(length(sufficiently_active_cells_indexes)) '/' num2str(length(active_bins)) ' sufficiently active cells, out of which ' num2str(length(significantly_tuned_and_active_cells_indexes)) ' cells (' num2str(round(100*length(significantly_tuned_and_active_cells_indexes)/length(sufficiently_active_cells_indexes))) '%) are significantly modulated by the encoded variable.'])
        else
            disp(['Found ' num2str(length(sufficiently_active_cells_indexes)) '/' num2str(length(active_bins)) ' sufficiently active cells'])
        end
        
        % Calculating the naive and shuffle information as a function of sample size:
        subsampling_repetitions=settings.subsampling_repetitions; % number of repetitions in the subsampling of the data
        T=size(spike_train,1); % total number of samples
        subsample_size=settings.subsample_fraction*T; % subsamples with size of different fractions of the data
        
        disp('Computing information as a function of subsample size:')
        if measures_to_estimate(1) || measures_to_estimate(2)
            if measures_to_estimate(3) % Compute SI and MI versus sample size
                [SI_naive_bit_spike_versus_sample_size,SI_shuffle_bit_spike_versus_sample_size,SI_naive_bit_sec_versus_sample_size,SI_shuffle_bit_sec_versus_sample_size,MI_naive_versus_sample_size,MI_shuffle_versus_sample_size]=...
                    compute_information_versus_sample_size(spike_train_significantly_tuned_cells,stimulus_trace,subsample_size,dt,subsampling_repetitions,measures_to_estimate);
            else % Compute only SI versus sample size
                [SI_naive_bit_spike_versus_sample_size,SI_shuffle_bit_spike_versus_sample_size,SI_naive_bit_sec_versus_sample_size,SI_shuffle_bit_sec_versus_sample_size]=...
                    compute_information_versus_sample_size(spike_train_significantly_tuned_cells,stimulus_trace,subsample_size,dt,subsampling_repetitions,measures_to_estimate);
            end
        elseif measures_to_estimate(3) % Compute only MI versus sample size
            [~,~,~,~,MI_naive_versus_sample_size,MI_shuffle_versus_sample_size]=...
                compute_information_versus_sample_size(spike_train_significantly_tuned_cells,stimulus_trace,subsample_size,dt,subsampling_repetitions,measures_to_estimate);
        end
        
        % Correcting the bias using the SSR BAE methods:
        plot_results=settings.plot_results;
        save_figures=settings.save_figures;
        figures_directory=settings.figures_directory;
        if measures_to_estimate(1) % Correcting the bias for SI in bit/spike
            
            units='bit/spike';
            
            % SSR method:
            [SI_SSR_bit_spike,average_SI_SSR_bit_spike,SI_SSR_stability_bit_spike,average_SI_SSR_stability_bit_spike]...
                =perform_SSR(SI_naive_bit_spike_versus_sample_size,SI_shuffle_bit_spike_versus_sample_size,subsample_size,units,plot_results,save_figures,figures_directory);
            
            % BAE method:
            [SI_BAE_bit_spike,average_SI_BAE_bit_spike,SI_BAE_fit_R_2_bit_spike,average_SI_BAE_fit_R_2_bit_spike]...
                =perform_BAE(SI_naive_bit_spike_versus_sample_size,subsample_size,units,plot_results,save_figures,figures_directory);
            
            SI_disagreement_bit_spike=SI_BAE_bit_spike-SI_SSR_bit_spike;
            average_SI_disagreement_bit_spike=average_SI_BAE_bit_spike-average_SI_SSR_bit_spike;
        end
        
        % Correcting the bias for SI in bit/sec:
        if measures_to_estimate(2)
            units='bit/sec';
            
            % SSR method:
            [SI_SSR_bit_sec,average_SI_SSR_bit_sec,SI_SSR_stability_bit_sec,average_SI_SSR_stability_bit_sec]...
                =perform_SSR(SI_naive_bit_sec_versus_sample_size,SI_shuffle_bit_sec_versus_sample_size,subsample_size,units,plot_results,save_figures,figures_directory);
            
            % BAE method:
            [SI_BAE_bit_sec,average_SI_BAE_bit_sec,SI_BAE_fit_R_2_bit_sec,average_SI_BAE_fit_R_2_bit_sec]...
                =perform_BAE(SI_naive_bit_sec_versus_sample_size,subsample_size,units,plot_results,save_figures,figures_directory);
            
            SI_disagreement_bit_sec=SI_BAE_bit_sec-SI_SSR_bit_sec;
            average_SI_disagreement_bit_sec=average_SI_BAE_bit_sec-average_SI_SSR_bit_sec;
        end
        
        % Correcting the bias for MI:
        if measures_to_estimate(3)
            units='bit';
            
            % SSR method:
            [MI_SSR,average_MI_SSR,MI_SSR_stability,average_MI_SSR_stability]...
                =perform_SSR(MI_naive_versus_sample_size,MI_shuffle_versus_sample_size,subsample_size,units,plot_results,save_figures,figures_directory);
            
            % BAE method:
            [MI_BAE,average_MI_BAE,MI_BAE_fit_R_2,average_MI_BAE_fit_R_2]...
                =perform_BAE(MI_naive_versus_sample_size,subsample_size,units,plot_results,save_figures,figures_directory);
            
            MI_disagreement=MI_BAE-MI_SSR;
            average_MI_disagreement=average_MI_BAE-average_MI_SSR;
        end
        
        % plotting the results:
        if plot_results || save_figures
            if measures_to_estimate(1) % for SI in bit/spike
                
                % Cross validation:
                if plot_results
                    figure
                else
                    figure('Visible','off')
                end
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
                if save_figures
                    savefig(fullfile(figures_directory,'BAE versus SSR - SI bit per spike.fig'))
                    saveas(gcf,fullfile(figures_directory,'BAE versus SSR - SI bit per spike'),'png')
                end
                
                % Estimation quality
                if plot_results
                    figure
                else
                    figure('Visible','off')
                end
                plot(SI_SSR_stability_bit_spike,SI_BAE_fit_R_2_bit_spike,'.','markersize',15,'color','b')
                xlim([0 1])
                ylim([0 1])
                axis square
                xlabel('SSR stability')
                ylabel('BAE fit R^2')
                title('SI (bit/spike)')
                set(gca,'fontsize',16)
                box off
                if save_figures
                    savefig(fullfile(figures_directory,'Estimation quality - SI bit per spike.fig'))
                    saveas(gcf,fullfile(figures_directory,'Estimation quality - SI bit per spike'),'png')
                end
            end
            
            if measures_to_estimate(2) % for SI in bit/sec
                
                % Cross validation:
                if plot_results
                    figure
                else
                    figure('Visible','off')
                end
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
                if save_figures
                    savefig(fullfile(figures_directory,'BAE versus SSR - SI bit per sec.fig'))
                    saveas(gcf,fullfile(figures_directory,'BAE versus SSR - SI bit per sec'),'png')
                end
                
                % Estimation quality
                if plot_results
                    figure
                else
                    figure('Visible','off')
                end
                plot(SI_SSR_stability_bit_sec,SI_BAE_fit_R_2_bit_sec,'.','markersize',15,'color','b')
                xlim([0 1])
                ylim([0 1])
                axis square
                xlabel('SSR stability')
                ylabel('BAE fit R^2')
                title('SI (bit/sec)')
                set(gca,'fontsize',16)
                box off
                if save_figures
                    savefig(fullfile(figures_directory,'Estimation quality - SI bit per sec.fig'))
                    saveas(gcf,fullfile(figures_directory,'Estimation quality - SI bit per sec'),'png')
                end
            end
            
            if measures_to_estimate(3) % for MI
                
                % Cross validation:
                if plot_results
                    figure
                else
                    figure('Visible','off')
                end
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
                if save_figures
                    savefig(fullfile(figures_directory,'BAE versus SSR - MI.fig'))
                    saveas(gcf,fullfile(figures_directory,'BAE versus SSR - MI'),'png')
                end
                
                % Estimation quality
                if plot_results
                    figure
                else
                    figure('Visible','off')
                end
                plot(MI_SSR_stability,MI_BAE_fit_R_2,'.','markersize',15,'color','b')
                xlim([0 1])
                ylim([0 1])
                axis square
                xlabel('SSR stability')
                ylabel('BAE fit R^2')
                title('MI')
                set(gca,'fontsize',16)
                box off
                if save_figures
                    savefig(fullfile(figures_directory,'Estimation quality - MI.fig'))
                    saveas(gcf,fullfile(figures_directory,'Estimation quality - MI'),'png')
                end
            end
        end
        
        % Saving the final results in a single data structure:
        % General parameters:
        unbiased_information_estimation_results.settings=struct;
        unbiased_information_estimation_results.settings.dt=dt;
        if estimate_only_significantly_tuned_cells
            unbiased_information_estimation_results.settings.num_shuffles=num_shuffles;
            unbiased_information_estimation_results.settings.shuffle_type=shuffle_type;
            unbiased_information_estimation_results.settings.tuning_significance_threshold=tuning_significance_threshold;
            unbiased_information_estimation_results.firing_statistics.p_value_significantly_tuned_and_active_cells=p_value_significantly_tuned_and_active_cells;
        end
        unbiased_information_estimation_results.settings.active_bins_threshold=active_bins_threshold;
        unbiased_information_estimation_results.settings.firing_rate_threshold=firing_rate_threshold;
        unbiased_information_estimation_results.settings.subsampling_repetitions=subsampling_repetitions;
        
        % Data statistics:
        unbiased_information_estimation_results.firing_statistics=struct;
        unbiased_information_estimation_results.firing_statistics.fraction_sufficiently_active_cells=fraction_sufficiently_active_cells;
        unbiased_information_estimation_results.firing_statistics.fraction_significantly_tuned_and_active_cells=fraction_significantly_tuned_and_active_cells;
        unbiased_information_estimation_results.firing_statistics.fraction_of_significantly_tuned_from_active_cells=fraction_of_significantly_tuned_from_active_cells;
        unbiased_information_estimation_results.firing_statistics.average_rates_significantly_tuned_cells=average_rates_significantly_tuned_cells;
        unbiased_information_estimation_results.firing_statistics.active_bins_significantly_tuned_cells=active_bins_significantly_tuned_cells;
        unbiased_information_estimation_results.firing_statistics.sufficiently_active_cells_indexes=sufficiently_active_cells_indexes;
        unbiased_information_estimation_results.firing_statistics.significantly_tuned_and_active_cells_indexes=significantly_tuned_and_active_cells_indexes;
        
        % Estimated information:
        if measures_to_estimate(1) % SI in bit/spike
            % For individual cells:
            unbiased_information_estimation_results.information.SI_naive_bit_spike=SI_naive_bit_spike_significantly_tuned_cells;
            unbiased_information_estimation_results.information.SI_SSR_bit_spike=SI_SSR_bit_spike;
            unbiased_information_estimation_results.information.SI_BAE_bit_spike=SI_BAE_bit_spike;
            unbiased_information_estimation_results.information.SI_disagreement_bit_spike=SI_disagreement_bit_spike;
            unbiased_information_estimation_results.information.SI_SSR_stability_bit_spike=SI_SSR_stability_bit_spike;
            unbiased_information_estimation_results.information.SI_BAE_fit_R_2_bit_spike=SI_BAE_fit_R_2_bit_spike;
            
            % Average across the population:
            unbiased_information_estimation_results.information.average_SI_naive_bit_spike=average_SI_naive_bit_spike_significantly_tuned_cells;
            unbiased_information_estimation_results.information.average_SI_SSR_bit_spike=average_SI_SSR_bit_spike;
            unbiased_information_estimation_results.information.average_SI_BAE_bit_spike=average_SI_BAE_bit_spike;
            unbiased_information_estimation_results.information.average_SI_disagreement_bit_spike=average_SI_disagreement_bit_spike;
            unbiased_information_estimation_results.information.average_SI_SSR_stability_bit_spike=average_SI_SSR_stability_bit_spike;
            unbiased_information_estimation_results.information.average_SI_BAE_fit_R_2_bit_spike=average_SI_BAE_fit_R_2_bit_spike;
            
            % Average across spikes (weighted by the cells firing rates):
            unbiased_information_estimation_results.information.weighted_average_SI_naive_bit_spike=sum(SI_naive_bit_spike_significantly_tuned_cells.*average_rates_significantly_tuned_cells')/sum(average_rates_significantly_tuned_cells);
            unbiased_information_estimation_results.information.weighted_average_SI_SSR_bit_spike=sum(SI_SSR_bit_spike.*average_rates_significantly_tuned_cells')/sum(average_rates_significantly_tuned_cells);
            unbiased_information_estimation_results.information.weighted_average_SI_BAE_bit_spike=sum(SI_BAE_bit_spike.*average_rates_significantly_tuned_cells')/sum(average_rates_significantly_tuned_cells);
            unbiased_information_estimation_results.information.weighted_average_SI_disagreement_bit_spike=sum(SI_disagreement_bit_spike.*average_rates_significantly_tuned_cells')/sum(average_rates_significantly_tuned_cells);
        end
        
        if measures_to_estimate(2) % SI in bit/sec
            % For individual cells:
            unbiased_information_estimation_results.information.SI_naive_bit_sec=SI_naive_bit_sec_significantly_tuned_cells;
            unbiased_information_estimation_results.information.SI_SSR_bit_sec=SI_SSR_bit_sec;
            unbiased_information_estimation_results.information.SI_BAE_bit_sec=SI_BAE_bit_sec;
            unbiased_information_estimation_results.information.SI_disagreement_bit_sec=SI_disagreement_bit_sec;
            unbiased_information_estimation_results.information.SI_SSR_stability_bit_sec=SI_SSR_stability_bit_sec;
            unbiased_information_estimation_results.information.SI_BAE_fit_R_2_bit_sec=SI_BAE_fit_R_2_bit_sec;
            
            % Average across the population:
            unbiased_information_estimation_results.information.average_SI_naive_bit_sec=average_SI_naive_bit_sec_significantly_tuned_cells;
            unbiased_information_estimation_results.information.average_SI_SSR_bit_sec=average_SI_SSR_bit_sec;
            unbiased_information_estimation_results.information.average_SI_BAE_bit_sec=average_SI_BAE_bit_sec;
            unbiased_information_estimation_results.information.average_SI_disagreement_bit_sec=average_SI_disagreement_bit_sec;
            unbiased_information_estimation_results.information.average_SI_SSR_stability_bit_sec=average_SI_SSR_stability_bit_sec;
            unbiased_information_estimation_results.information.average_SI_BAE_fit_R_2_bit_sec=average_SI_BAE_fit_R_2_bit_sec;
        end
        
        if measures_to_estimate(3) % MI
            % For individual cells:
            unbiased_information_estimation_results.information.MI_naive=MI_naive_significantly_tuned_cells;
            unbiased_information_estimation_results.information.MI_SSR=MI_SSR;
            unbiased_information_estimation_results.information.MI_BAE=MI_BAE;
            unbiased_information_estimation_results.information.MI_disagreement=MI_disagreement;
            unbiased_information_estimation_results.information.MI_SSR_stability=MI_SSR_stability;
            unbiased_information_estimation_results.information.MI_BAE_fit_R_2=MI_BAE_fit_R_2;
            
            % Average across the population:
            unbiased_information_estimation_results.information.average_MI_naive=average_MI_naive_significantly_tuned_cells;
            unbiased_information_estimation_results.information.average_MI_SSR=average_MI_SSR;
            unbiased_information_estimation_results.information.average_MI_BAE=average_MI_BAE;
            unbiased_information_estimation_results.information.average_MI_disagreement=average_MI_disagreement;
            unbiased_information_estimation_results.information.average_MI_SSR_stability=average_MI_SSR_stability;
            unbiased_information_estimation_results.information.average_MI_BAE_fit_R_2=average_MI_BAE_fit_R_2;
        end
        disp('Finished analyzing data set')
        
    else
        disp(['Found ' num2str(length(sufficiently_active_cells_indexes)) '/' num2str(length(active_bins)) ' sufficiently active cells'])
        disp('No significantly tuned cells were found')
    end
else
    disp('No sufficiently active cells were found')
end

end
