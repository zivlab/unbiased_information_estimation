function [SSR_information,mean_SSR_information,SSR_stability,mean_SSR_stability]=perform_SSR_simulated_data(information_versus_sample_size,shuffle_information_versus_sample_size,true_information,subsample_size,units,figures_directory,plot_results)
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
% 3. true_information - Vector of size N with the true information for each neuron.
% 4. subsample_size - Vector of T different sample sizes
% 5. units - Either bit/spike, bit/sec, or bit
% 6. figures_directory - Path to save the results figure
% 7. plot_results - 1 for plotting and 0 for not

% Outputs:
% 1. SSR_information - Vector of size N with the estimated
% information for each neuron.
% 2. mean_SSR_information - Single value with the estimated
% average information for the mean..

% computing the SSR information:
mean_information_versus_sample_size=mean(information_versus_sample_size,'omitnan');
mean_shuffle_information_versus_sample_size=mean(shuffle_information_versus_sample_size,'omitnan');
SSR_information_versus_sample_size=information_versus_sample_size-shuffle_information_versus_sample_size.*(information_versus_sample_size(:,1)-information_versus_sample_size)./(shuffle_information_versus_sample_size(:,1)-shuffle_information_versus_sample_size);
SSR_information=information_versus_sample_size(:,end)-shuffle_information_versus_sample_size(:,end).*(information_versus_sample_size(:,1)-information_versus_sample_size(:,end))./(shuffle_information_versus_sample_size(:,1)-shuffle_information_versus_sample_size(:,end));
mean_information_short_duration=mean(information_versus_sample_size(:,1),'omitnan'); % for subsample duration t1
mean_shuffle_information_short_duration=mean(shuffle_information_versus_sample_size(:,1),'omitnan');
mean_SSR_information_versus_sample_size=mean_information_versus_sample_size-mean_shuffle_information_versus_sample_size.*(mean_information_short_duration-mean_information_versus_sample_size)./(mean_shuffle_information_short_duration-mean_shuffle_information_versus_sample_size);
mean_SR_information_versus_sample_size=mean_information_versus_sample_size-mean_shuffle_information_versus_sample_size;
mean_SSR_information=mean_SSR_information_versus_sample_size(end);
SSR_stability=1-std(SSR_information_versus_sample_size,0,2,'omitnan')./SSR_information.*subsample_size(2)./subsample_size(end);
mean_SSR_stability=1-std(mean_SSR_information_versus_sample_size,'omitnan')./mean_SSR_information.*subsample_size(2)./subsample_size(end);

% plotting the mean average results for the SSR method:
if plot_results
    figure;
    plot((1:length(mean_information_versus_sample_size))./length(mean_information_versus_sample_size),mean_information_versus_sample_size,'color','b','linewidth',2)
    hold on
    plot((1:length(mean_information_versus_sample_size))./length(mean_information_versus_sample_size),mean_shuffle_information_versus_sample_size,'color','k','linewidth',2)
    plot((1:length(mean_information_versus_sample_size))./length(mean_information_versus_sample_size),mean_SR_information_versus_sample_size,'-r','linewidth',2)
    plot((1:length(mean_information_versus_sample_size))./length(mean_information_versus_sample_size),mean_SSR_information_versus_sample_size,'-m','linewidth',2)
    plot([0 1],[mean(true_information,'omitnan') mean(true_information,'omitnan')],'-c','linewidth',2)
    plot((1:length(mean_information_versus_sample_size))./length(mean_information_versus_sample_size),mean_information_versus_sample_size,'color','b','linewidth',2)
    plot((1:length(mean_information_versus_sample_size))./length(mean_information_versus_sample_size),mean_shuffle_information_versus_sample_size,'color','k','linewidth',2)
    plot((1:length(mean_information_versus_sample_size))./length(mean_information_versus_sample_size),mean_SR_information_versus_sample_size,'-r','linewidth',2)
    plot((1:length(mean_information_versus_sample_size))./length(mean_information_versus_sample_size),mean_SSR_information_versus_sample_size,'-m','linewidth',2)
    ylim([0 1.1*max(mean_information_versus_sample_size)])
    xlim([0 1])
    if strcmp(units,'bit')
        ylim([0 1.1*max(mean_information_versus_sample_size)])
    end
    xlabel('Subsample fraction')
    if strcmp(units,'bit/spike') || strcmp(units,'bit/sec')
        ylabel(['SI (' units ')'])
    else
        ylabel(['MI (' units ')'])
    end
    legend('Naive','Shuffle','SR','SSR','True')
    legend boxoff
    set(gca,'fontsize',16)
    box off
    axis square
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
    
    % plotting the individual cells results for the SSR method:
    figure
    plot(true_information,SSR_information,'.','markersize',15,'color','m')
    hold on
    plot([0 1.1*max(information_versus_sample_size(:,end))],[0 1.1*max(information_versus_sample_size(:,end))],'--k','linewidth',2)
    xlim([0 1.1*max(information_versus_sample_size(:,end))])
    ylim([0 1.1*max(information_versus_sample_size(:,end))])
    axis square
    if strcmp(units,'bit/spike') || strcmp(units,'bit/sec')
        xlabel(['True SI (' units ')'])
    else
        xlabel(['True MI (' units ')'])
    end
    ylabel(['SSR estimation (' units ')'])
    set(gca,'fontsize',16)
    box off
    if strcmp(units,'bit/spike')
        savefig(fullfile(figures_directory,'SSR versus true information - SI bit per spike.fig'))
        saveas(gcf,fullfile(figures_directory,'SSR versus true information - SI bit per spike'),'png')
    elseif strcmp(units,'bit/sec')
        savefig(fullfile(figures_directory,'SSR versus true information - SI bit per sec.fig'))
        saveas(gcf,fullfile(figures_directory,'SSR versus true information - SI bit per sec'),'png')
    else
        savefig(fullfile(figures_directory,'SSR versus true information - MI.fig'))
        saveas(gcf,fullfile(figures_directory,'SSR versus true information - MI'),'png')
    end
end
