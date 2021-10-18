function [SSR_information,mean_SSR_information,SSR_stability,mean_SSR_stability]=perform_SSR(information_versus_sample_size,shuffle_information_versus_sample_size,subsample_size,units,stability_threshold,mean_stability_threshold,figures_directory)
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
% 5. stability_threshold - warning for unstable estiamtion for individual cells 
% 6. mean_stability_threshold - warning for unstable estiamtion for the mean
% 7. figures_directory - Path to save the results figure

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
figure; 
plot((1:length(mean_information_versus_sample_size))./length(mean_information_versus_sample_size),mean_information_versus_sample_size,'color','b','linewidth',2)
hold on
plot((1:length(mean_information_versus_sample_size))./length(mean_information_versus_sample_size),mean_shuffle_information_versus_sample_size,'color','k','linewidth',2)
plot((1:length(mean_information_versus_sample_size))./length(mean_information_versus_sample_size),mean_SR_information_versus_sample_size,'-r','linewidth',2)
plot((1:length(mean_information_versus_sample_size))./length(mean_information_versus_sample_size),mean_SSR_information_versus_sample_size,'-m','linewidth',2)
ylim([0 ceil(max(mean_information_versus_sample_size))])
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
legend('Naive','Shuffle','SR','SSR')
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

% checking if mean bias correction is sufficiently stable:
if mean_SSR_stability<mean_stability_threshold
    warndlg(['The SSR method is not very stable (mean stability<' num2str(mean_stability_threshold) '). Try choosing a larger minimal subsample size.'])
end

% plotting the individual cells results for the SSR method:
figure
plot(information_versus_sample_size(:,end),SSR_information,'.','markersize',15,'color','m')
hold on
plot(information_versus_sample_size(SSR_stability<stability_threshold,end),SSR_information(SSR_stability<stability_threshold),'.','markersize',15,'color',[0.5 0.5 0.5])
plot([0 1.1*max(information_versus_sample_size(:,end))],[0 1.1*max(information_versus_sample_size(:,end))],'--k','linewidth',2)
xlim([0 1.1*max(information_versus_sample_size(:,end))])
ylim([0 1.1*max(information_versus_sample_size(:,end))])
axis square
if strcmp(units,'bit/spike') || strcmp(units,'bit/sec')
    xlabel(['Naive SI (' units ')'])
    ylabel(['SSR SI (' units ')'])
else
    xlabel(['Naive MI (' units ')'])
    ylabel(['SSR MI (' units ')'])
end
ylabel(['SSR (' units ')'])
set(gca,'fontsize',16)
box off
if sum(SSR_stability<stability_threshold)>0
    legend('Stable estimation','Unstable estimation','Location','northwest')
    legend boxoff
end
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

% finding cells with unstable bias corrections:
if sum(SSR_stability<stability_threshold)>0
    warndlg(['The SSR method is less stable in ' num2str(sum(SSR_stability<stability_threshold)) '/' num2str(length(SSR_stability)) ' cells (stability<' num2str(stability_threshold) '). Consider discarding these cells from the analysis.'])
end
