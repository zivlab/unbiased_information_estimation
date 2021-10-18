function [BAE_information,mean_BAE_information,BAE_fit_R_2,BAE_fit_R_2_mean]=perform_BAE_simulated_data(information_versus_sample_size,true_information,subsample_size,units,BAE_fit_R_2_threshold,mean_BAE_fit_R_2_threshold,figures_directory)
% This functions corrects the upward bias in the naive calculation of
% information content for limited sample sizes using the bounded extrapolation (BAE) method.
% BAE is based on fitting the function of how the information changes with sample size and extrapolating it to infinity.
% The obtained results are compared against the unbounded extrapolation (AE) method. 

% Inputs:
% 1. information_versus_sample_size - Matrix of size TxN with the estimated
% information of each of N neurons as a function of T different sample
% sizes.
% 2. true_information - Vector of size N with the true information for each neuron.
% 3. subsample_size - Vector of T different sample sizes
% 4. units - Either bit/spike, bit/sec, or bit
% 5. BAE_fit_R_2_threshold - warning for unstable estiamtion for individual cells 
% 6. populaiton_BAE_fit_R_2_threshold - warning for unstable estiamtion for the mean
% 7. figures_directory - Path to save the results figure

% Outputs:
% 1. BAE_information - Vector of size N with the estimated
% information for each neuron.
% 2. mean_BAE_information - Single value with the estimated
% average information for the mean.

% extrapolating the information for the mean average with information(t)=a+b/(1+ct) - BAE:
mean_information_versus_sample_size=mean(information_versus_sample_size,'omitnan');
a_0=mean(information_versus_sample_size(:,end),'omitnan');
b_0=mean(information_versus_sample_size(:,1),'omitnan')-mean(information_versus_sample_size(:,end),'omitnan');
middle_index=round(length(subsample_size)/2);
middle_sample_size=subsample_size(middle_index);
c_0=1./middle_sample_size*(mean(information_versus_sample_size(:,1),'omitnan')-mean(information_versus_sample_size(:,middle_index),'omitnan'))/(mean(information_versus_sample_size(:,middle_index),'omitnan')-mean(information_versus_sample_size(:,end),'omitnan'));
initial_parameters=[a_0 b_0 c_0];
F_model = @(x,xdata)...
    x(1)+(x(2))./(1+x(3).*xdata);
lb = [0 0 0];
ub = [Inf Inf Inf];
options = statset('MaxIter',1000, 'MaxFunEvals',2000);

% finding the parameters that best fit the data:
mean_BAE_fit_params=lsqcurvefit(F_model,initial_parameters,subsample_size,mean_information_versus_sample_size,lb,ub,options);
BAE_fitted_model=mean_BAE_fit_params(1)+(mean_BAE_fit_params(2))./(1+mean_BAE_fit_params(3).*subsample_size);
mean_BAE_information=mean_BAE_fit_params(1);
BAE_fit_R_2_mean=1-mean((BAE_fitted_model-mean_information_versus_sample_size).^2)./var(mean_information_versus_sample_size);

% extrapolating the information for the mean average with information(t)=a+b/t+c/t^2 - AE:
a_0=mean(information_versus_sample_size(:,end),'omitnan');
b_0=mean(information_versus_sample_size(:,1),'omitnan')-mean(information_versus_sample_size(:,end),'omitnan');
c_0=mean(information_versus_sample_size(:,1),'omitnan')-mean(information_versus_sample_size(:,end),'omitnan');
initial_parameters=[a_0 b_0 c_0];
F_model = @(x,xdata)...
    x(1)+x(2)./xdata+x(3)./xdata.^2;
lb = [0 0 0];
ub = [Inf Inf inf];
options = statset('MaxIter',1000, 'MaxFunEvals',2000);

% finding the parameters that best fit the data:
AE_fit_params=lsqcurvefit(F_model,initial_parameters,subsample_size,mean_information_versus_sample_size,lb,ub,options);
AE_fitted_model=AE_fit_params(1)+AE_fit_params(2)./subsample_size+AE_fit_params(3)./subsample_size.^2;

% extrapolating the information for each cell:
N=size(information_versus_sample_size,1);
BAE_information=nan(N,1);
BAE_fit_R_2=nan(N,1);
for n=1:N
    this_information_versus_sample_size=information_versus_sample_size(n,:);
    a_0=this_information_versus_sample_size(end);
    b_0=this_information_versus_sample_size(1)-this_information_versus_sample_size(end);
    c_0=1./middle_sample_size*(this_information_versus_sample_size(1)-this_information_versus_sample_size(middle_index))/(this_information_versus_sample_size(middle_index)-this_information_versus_sample_size(end));
    
    initial_parameters=[a_0 b_0 c_0];
    F_model = @(x,xdata)...
        x(1)+x(2)./(1+x(3).*xdata);
    lb = [0 0 0];
    ub = [Inf Inf Inf];
    options = statset('MaxIter',1000, 'MaxFunEvals',2000);
    % finding the parameters that best fit the data:
    this_BAE_fit_params=lsqcurvefit(F_model,initial_parameters,subsample_size,this_information_versus_sample_size,lb,ub,options);
    BAE_information(n)=this_BAE_fit_params(1);
    this_BAE_fitted_model=this_BAE_fit_params(1)+(this_BAE_fit_params(2))./(1+this_BAE_fit_params(3).*subsample_size);
    BAE_fit_R_2(n)=1-mean((this_BAE_fitted_model-this_information_versus_sample_size).^2)./var(this_information_versus_sample_size);
end

% plotting the mean average results for the extrapolation method:
figure
plot(subsample_size/subsample_size(end),mean_information_versus_sample_size,'ob','linewidth',2)
hold on
plot(subsample_size/subsample_size(end),AE_fitted_model,'-','color',[1 0.5 0],'linewidth',2)
plot([0 1],[AE_fit_params(1) AE_fit_params(1)],'--','color',[1 0.5 0],'linewidth',2)
plot(subsample_size/subsample_size(end),BAE_fitted_model,'-g','linewidth',2)
plot([0 1],[mean_BAE_fit_params(1) mean_BAE_fit_params(1)],'--g','linewidth',2)
plot([0 1],[mean(true_information,'omitnan') mean(true_information,'omitnan')],'-c','linewidth',2)
plot(subsample_size/subsample_size(end),AE_fitted_model,'-','color',[1 0.5 0],'linewidth',2)
plot([0 1],[AE_fit_params(1) AE_fit_params(1)],'--','color',[1 0.5 0],'linewidth',2)
plot(subsample_size/subsample_size(end),BAE_fitted_model,'-g','linewidth',2)
plot([0 1],[mean_BAE_fit_params(1) mean_BAE_fit_params(1)],'--g','linewidth',2)
plot(subsample_size/subsample_size(end),mean_information_versus_sample_size,'ob','linewidth',2)
xlim([0 1])
ylim([0 1.1*(mean_information_versus_sample_size(1))])
if strcmp(units,'bit')
    ylim([0 1.1*mean_information_versus_sample_size(1)])
end
xlabel('Subsample fraction')
if strcmp(units,'bit/spike') || strcmp(units,'bit/sec')
    ylabel(['SI (' units ')'])
else
    ylabel(['MI (' units ')'])
end
legend('Naive','AE fitted model','AE estimation','BAE fitted model','BAE estimation','True')
legend boxoff
set(gca,'fontsize',16)
box off
axis square
if strcmp(units,'bit/spike')
savefig(fullfile(figures_directory,'BAE method - SI bit per spike.fig'))
saveas(gcf,fullfile(figures_directory,'BAE method - SI bit per spike'),'png')
elseif strcmp(units,'bit/sec')
    savefig(fullfile(figures_directory,'BAE method - SI bit per sec.fig'))
    saveas(gcf,fullfile(figures_directory,'BAE method - SI bit per sec'),'png')
else
    savefig(fullfile(figures_directory,'BAE method - MI.fig'))
    saveas(gcf,fullfile(figures_directory,'BAE method - MI'),'png')
end

% checking if mean fit accuracy is sufficiently high:
if BAE_fit_R_2_mean<mean_BAE_fit_R_2_threshold
    warndlg(['The fit of the BAE method is not very accurate (mean R^2<' num2str(mean_BAE_fit_R_2_threshold) '). Try choosing a larger minimal subsample size.'])
end

% plotting the cell-level results for the extrapolation method:
figure
plot(true_information,BAE_information,'.','markersize',15,'color','g')
hold on
plot(true_information(BAE_fit_R_2<BAE_fit_R_2_threshold),BAE_information(BAE_fit_R_2<BAE_fit_R_2_threshold),'.','markersize',15,'color',[0.5 0.5 0.5])
plot([0 1.1*max(information_versus_sample_size(:,end))],[0 1.1*max(information_versus_sample_size(:,end))],'--k','linewidth',2)
xlim([0 1.1*max(information_versus_sample_size(:,end))])
ylim([0 1.1*max(information_versus_sample_size(:,end))])
axis square
if strcmp(units,'bit/spike') || strcmp(units,'bit/sec')
    xlabel(['Naive SI (' units ')'])
    ylabel(['BAE SI (' units ')'])
else
    xlabel(['Naive MI (' units ')'])
    ylabel(['BAE MI(' units ')'])
end
ylabel(['BAE (' units ')'])
set(gca,'fontsize',16)
box off
if sum(BAE_fit_R_2<BAE_fit_R_2_threshold)>0
    legend('Accurate fit','Inaccurate fit','Location','northwest')
    legend boxoff
end
if strcmp(units,'bit/spike')
    savefig(fullfile(figures_directory,'BAE versus true information - SI bit per spike.fig'))
    saveas(gcf,fullfile(figures_directory,'BAE versus true information - SI bit per spike'),'png')
elseif strcmp(units,'bit/sec')
    savefig(fullfile(figures_directory,'BAE versus true information - SI bit per sec.fig'))
    saveas(gcf,fullfile(figures_directory,'BAE versus true information - SI bit per sec'),'png')
else
    savefig(fullfile(figures_directory,'BAE versus true information - MI.fig'))
    saveas(gcf,fullfile(figures_directory,'BAE versus true information - MI'),'png')
end

% finding cells with inaccurate fit:
if sum(BAE_fit_R_2<BAE_fit_R_2_threshold)>0
    warndlg(['The fit of the BAE method is less accurate in ' num2str(sum(BAE_fit_R_2<BAE_fit_R_2_threshold)) '/' num2str(N) ' cells (R^2<' num2str(BAE_fit_R_2_threshold) '). Consider discarding these cells from the analysis.'])
end