function [SSR_information,average_SSR_information,SSR_stability,average_SSR_stability]=perform_SSR(information_versus_sample_size,shuffle_information_versus_sample_size,subsample_size,units,plot_results,save_figures,figures_directory)
% This functions corrects the upward bias in the naive calculation of
% information for limited sample sizes using the scaled shuffle reduction (SSR) method. 
% SSR is based on assuming a fixed bias ratio between the naive and shuffle information 
% and using two different subsample sizes to find this ratio and subtract a the scaled shuffle information.
% The obtained results are compared against the shuffle reduction (SR) method. 

% Inputs:
% 1. information_versus_sample_size - Matrix of size TxN with the estimated
% information of each of N neurons as a function of T different sample sizes
% 2. shuffle_information_versus_sample_size - Matrix of size TxN with the estimated
% bounds and 3 when using only the data
% 3. subsample_size - Vector of T different sample sizes
% 4. units - Either bit/spike, bit/sec, or bit
% 5. plot_results - 1 for plotting and 0 for not
% 6. save_figures - 1 for saving the figures and 0 for not
% 7. figures_directory - Path for saving the figures

% Outputs:
% 1. SSR_information - Vector of size N with the estimated
% information for each neuron.
% 2. average_SSR_information - Single value with the estimated
% average information for the population.
% 3. SSR_stability - Vector of size N with the SSR stability for each neuron.
% 4. average_SSR_stability - Single value with the average SSR stability for the population.

% computing the SSR information:

SSR_information_versus_sample_size=information_versus_sample_size-shuffle_information_versus_sample_size.*(information_versus_sample_size(:,1)-information_versus_sample_size)./(shuffle_information_versus_sample_size(:,1)-shuffle_information_versus_sample_size);
SSR_information=information_versus_sample_size(:,end)-shuffle_information_versus_sample_size(:,end).*(information_versus_sample_size(:,1)-information_versus_sample_size(:,end))./(shuffle_information_versus_sample_size(:,1)-shuffle_information_versus_sample_size(:,end));
SSR_stability=1-std(SSR_information_versus_sample_size,0,2,'omitnan')./SSR_information.*subsample_size(2)./subsample_size(end);

if size(information_versus_sample_size,1)>1
    average_information_versus_sample_size=mean(information_versus_sample_size,'omitnan');
    average_shuffle_information_versus_sample_size=mean(shuffle_information_versus_sample_size,'omitnan');
    average_information_short_duration=mean(information_versus_sample_size(:,1),'omitnan'); % for subsample duration t1
    average_shuffle_information_short_duration=mean(shuffle_information_versus_sample_size(:,1),'omitnan');
    average_SSR_information_versus_sample_size=average_information_versus_sample_size-average_shuffle_information_versus_sample_size.*(average_information_short_duration-average_information_versus_sample_size)./(average_shuffle_information_short_duration-average_shuffle_information_versus_sample_size);
    average_SR_information_versus_sample_size=average_information_versus_sample_size-average_shuffle_information_versus_sample_size;
    average_SSR_information=average_SSR_information_versus_sample_size(end);
    average_SSR_stability=1-std(average_SSR_information_versus_sample_size,'omitnan')./average_SSR_information.*subsample_size(2)./subsample_size(end);
else
    average_information_versus_sample_size=information_versus_sample_size;
    average_shuffle_information_versus_sample_size=shuffle_information_versus_sample_size;
    average_SR_information_versus_sample_size=average_information_versus_sample_size-average_shuffle_information_versus_sample_size;
    average_SSR_information_versus_sample_size=SSR_information_versus_sample_size;
    average_SSR_information=SSR_information;
    average_SSR_stability=SSR_stability;
end

% plotting the average average results for the SSR method:
if plot_results || save_figures
    if plot_results
        figure
    else
        figure('Visible','off')
    end
    plot((1:length(average_information_versus_sample_size))./length(average_information_versus_sample_size),average_information_versus_sample_size,'color','b','linewidth',2)
    hold on
    plot((1:length(average_information_versus_sample_size))./length(average_information_versus_sample_size),average_shuffle_information_versus_sample_size,'color','k','linewidth',2)
    plot((1:length(average_information_versus_sample_size))./length(average_information_versus_sample_size),average_SR_information_versus_sample_size,'-r','linewidth',2)
    plot((1:length(average_information_versus_sample_size))./length(average_information_versus_sample_size),average_SSR_information_versus_sample_size,'-m','linewidth',2)
    ylim([0 ceil(max(average_information_versus_sample_size))])
    xlim([0 1])
    if strcmp(units,'bit')
        ylim([0 1.1*max(average_information_versus_sample_size)])
    end
    xlabel('Subsample fraction')
    if strcmp(units,'bit/spike') || strcmp(units,'bit/sec')
        ylabel(['SI (' units ')'])
    else
        ylabel(['MI (' units ')'])
    end
    legend('Naive','Shuffle','SR','SSR')
    legend boxoff
    set(gca,'fontsize',16)
    box off
    axis square
    if save_figures
        if strcmp(units,'bit/spike')
            savefig(fullfile(figures_directory,'SSR method - SI bit per spike.fig'))
            saveas(gcf,fullfile(figures_directory,'SSR method - SI bit per spike'),'png')
        elseif strcmp(units,'bit/sec')
            savefig(fullfile(figures_directory,'SSR method - SI bit per sec.fig'))
            saveas(gcf,fullfile(figures_directory,'SSR method - SI bit per sec'),'png')
        else
            savefig(fullfile(figures_directory,'SSR method - MI.fig'))
            saveas(gcf,fullfile(figures_directory,'SSR method - MI'),'png')
        end
    end
    % plotting the individual cells results for the SSR method:
    if plot_results
        figure
    else
        figure('Visible','off')
    end
    plot(information_versus_sample_size(:,end),SSR_information,'.','markersize',15,'color','m')
    hold on
    plot([0 1.1*max(information_versus_sample_size(:,end))],[0 1.1*max(information_versus_sample_size(:,end))],'--k','linewidth',2)
    xlim([0 1.1*max(information_versus_sample_size(:,end))])
    ylim([0 1.1*max(information_versus_sample_size(:,end))])
    axis square
    if strcmp(units,'bit/spike') || strcmp(units,'bit/sec')
        xlabel(['Naive SI (' units ')'])
    else
        xlabel(['Naive MI (' units ')'])
    end
    ylabel(['SSR estimation (' units ')'])
    set(gca,'fontsize',16)
    box off
    if save_figures
        if strcmp(units,'bit/spike')
            savefig(fullfile(figures_directory,'SSR versus naive information - SI bit per spike.fig'))
            saveas(gcf,fullfile(figures_directory,'SSR versus naive information - SI bit per spike'),'png')
        elseif strcmp(units,'bit/sec')
            savefig(fullfile(figures_directory,'SSR versus naive information - SI bit per sec.fig'))
            saveas(gcf,fullfile(figures_directory,'SSR versus naive information - SI bit per sec'),'png')
        else
            savefig(fullfile(figures_directory,'SSR versus naive information - MI.fig'))
            saveas(gcf,fullfile(figures_directory,'SSR versus naive information - MI'),'png')
        end
    end
end
