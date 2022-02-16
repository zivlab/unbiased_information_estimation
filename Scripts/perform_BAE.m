function [BAE_information,average_BAE_information,BAE_fit_R_2,average_BAE_fit_R_2]=perform_BAE(information_versus_sample_size,subsample_size,units,plot_results,save_figures,figures_directory)
% This functions corrects the upward bias in the naive calculation of
% information content for limited sample sizes using the bounded asymptotic extrapolation (BAE) method.
% BAE is based on fitting the function of how the information changes with sample size and extrapolating it to infinity.
% The obtained results are compared against the unbounded aymptotic extrapolation (AE) method.

% Inputs:
% 1. information_versus_sample_size - Matrix of size TxN with the estimated
% information of each of N neurons as a function of T different sample
% sizes.
% 2. subsample_size - Vector of T different sample sizes
% 3. units - Either bit/spike, bit/sec, or bit% 5. plot_results - 1 for plotting and 0 for not
% 4. plot_results - 1 for plotting and 0 for not
% 5. save_figures - 1 for saving the figures and 0 for not
% 6. figures_directory - Path for saving the figures


% Outputs:
% 1. BAE_information - Vector of size N with the estimated
% information for each neuron.
% 2. average_BAE_information - Single value with the estimated
% average information for the population.
% 3. BAE_fit_R_2 - Vector of size N with the squared residuals of the
% fitted curve for each neuron
% 4. average_BAE_fit_R_2 - Single value with the average squared residuals for the population.

% extrapolating the information for the average average with information(t)=a+b/(1+ct) - BAE:
middle_index=round(length(subsample_size)/2);
middle_sample_size=subsample_size(middle_index);
if size(information_versus_sample_size,1)>1
    average_information_versus_sample_size=mean(information_versus_sample_size,'omitnan');
    a_0=mean(information_versus_sample_size(:,end),'omitnan');
    b_0=mean(information_versus_sample_size(:,1),'omitnan')-mean(information_versus_sample_size(:,end),'omitnan');
    c_0=1./middle_sample_size*(mean(information_versus_sample_size(:,1),'omitnan')-mean(information_versus_sample_size(:,middle_index),'omitnan'))/(mean(information_versus_sample_size(:,middle_index),'omitnan')-mean(information_versus_sample_size(:,end),'omitnan'));
else
    average_information_versus_sample_size=information_versus_sample_size;
    a_0=information_versus_sample_size(1,end);
    b_0=information_versus_sample_size(1,1)-information_versus_sample_size(1,end);
    c_0=1./middle_sample_size*(information_versus_sample_size(1,1)-information_versus_sample_size(1,middle_index))/(information_versus_sample_size(1,middle_index)-information_versus_sample_size(1,end));
end
initial_parameters=[a_0 b_0 c_0];
F_model = @(x,xdata)...
    x(1)+(x(2))./(1+x(3).*xdata);
lb = [0 0 0];
ub = [Inf Inf Inf];
options = statset('MaxIter',1000, 'MaxFunEvals',2000);

% finding the parameters that best fit the data (BAE):
average_BAE_fit_params=lsqcurvefit(F_model,initial_parameters,subsample_size,average_information_versus_sample_size,lb,ub,options);
BAE_fitted_model=average_BAE_fit_params(1)+(average_BAE_fit_params(2))./(1+average_BAE_fit_params(3).*subsample_size);
average_BAE_information=average_BAE_fit_params(1);
average_BAE_fit_R_2=1-mean((BAE_fitted_model-average_information_versus_sample_size).^2)./var(average_information_versus_sample_size);

% extrapolating the information for the average average with information(t)=a+b/t+c/t^2 - AE:
if size(information_versus_sample_size,1)>1
    a_0=mean(information_versus_sample_size(:,end),'omitnan');
    b_0=mean(information_versus_sample_size(:,1),'omitnan')-mean(information_versus_sample_size(:,end),'omitnan');
    c_0=mean(information_versus_sample_size(:,1),'omitnan')-mean(information_versus_sample_size(:,end),'omitnan');
else
    a_0=information_versus_sample_size(1,end);
    b_0=information_versus_sample_size(1,1)-information_versus_sample_size(1,end);
    c_0=information_versus_sample_size(1,1)-information_versus_sample_size(1,end);
end
initial_parameters=[a_0 b_0 c_0];
F_model = @(x,xdata)...
    x(1)+x(2)./xdata+x(3)./xdata.^2;
lb = [0 0 0];
ub = [Inf Inf inf];
options = statset('MaxIter',1000, 'MaxFunEvals',2000);

% finding the parameters that best fit the data (AE):
AE_fit_params=lsqcurvefit(F_model,initial_parameters,subsample_size,average_information_versus_sample_size,lb,ub,options);
AE_fitted_model=AE_fit_params(1)+AE_fit_params(2)./subsample_size+AE_fit_params(3)./subsample_size.^2;

% extrapolating the information for each cell:
if size(information_versus_sample_size,1)>1
    N=size(information_versus_sample_size,1);
    BAE_information=nan(N,1);
    BAE_fit_R_2=nan(N,1);
    for n=1:N
        this_information_versus_sample_size=information_versus_sample_size(n,:);
        if max(this_information_versus_sample_size)>0
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
    end
else
    BAE_information=average_BAE_information;
    BAE_fit_R_2=average_BAE_fit_R_2;
end

% plotting the average average results for the extrapolation method:
if plot_results || save_figures
    if plot_results
        figure
    else
        figure('Visible','off')
    end
    plot(subsample_size/subsample_size(end),average_information_versus_sample_size,'ob','linewidth',2)
    hold on
    plot(subsample_size/subsample_size(end),AE_fitted_model,'-','color',[1 0.5 0],'linewidth',2)
    plot([0 1],[AE_fit_params(1) AE_fit_params(1)],'--','color',[1 0.5 0],'linewidth',2)
    plot(subsample_size/subsample_size(end),BAE_fitted_model,'-g','linewidth',2)
    plot([0 1],[average_BAE_fit_params(1) average_BAE_fit_params(1)],'--g','linewidth',2)
    plot(subsample_size/subsample_size(end),average_information_versus_sample_size,'ob','linewidth',2)
    xlim([0 1])
    ylim([0 ceil(average_information_versus_sample_size(1))])
    if strcmp(units,'bit')
        ylim([0 1.1*average_information_versus_sample_size(1)])
    end
    xlabel('Subsample fraction')
    if strcmp(units,'bit/spike') || strcmp(units,'bit/sec')
        ylabel(['SI (' units ')'])
    else
        ylabel(['MI (' units ')'])
    end
    legend('Naive','AE fitted model','AE estimation','BAE fitted model','BAE estimation')
    legend boxoff
    set(gca,'fontsize',16)
    box off
    axis square
    if save_figures
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
    end
    
    % plotting the cell-level results for the extrapolation method:
    if plot_results
        figure
    else
        figure('Visible','off')
    end
    plot(information_versus_sample_size(:,end),BAE_information,'.','markersize',15,'color','g')
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
    ylabel(['BAE estimation (' units ')'])
    set(gca,'fontsize',16)
    box off
    if save_figures
        if strcmp(units,'bit/spike')
            savefig(fullfile(figures_directory,'BAE versus naive information - SI bit per spike.fig'))
            saveas(gcf,fullfile(figures_directory,'BAE versus naive information - SI bit per spike'),'png')
        elseif strcmp(units,'bit/sec')
            savefig(fullfile(figures_directory,'BAE versus naive information - SI bit per sec.fig'))
            saveas(gcf,fullfile(figures_directory,'BAE versus naive information - SI bit per sec'),'png')
        else
            savefig(fullfile(figures_directory,'BAE versus naive information - MI.fig'))
            saveas(gcf,fullfile(figures_directory,'BAE versus naive information - MI'),'png')
        end
    end
end
