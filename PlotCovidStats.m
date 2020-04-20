%% Setup
close all;
clearvars;

MIN_INF_THRESH = 1000;
FINAL_PLOT_P = 0.99;
DATE_FORMAT = 30;
LIN_SEARCH_PTS_PER_DIM = 100;
LIN_SEARCH_SCALE = 5;
MIN_POP_FOR_JOINED_PLOTS = 30e6;
PLOT_DIR_NAME = "Plots/" + datestr(now,DATE_FORMAT) + "/";
LATEST_DIR_NAME = "Latest/";

CountriesToTrack = fileread("CountriesToTrack.txt");
CountriesToTrack = split(CountriesToTrack,newline);
CountriesToTrack = string(CountriesToTrack(1:end-1));



%% Prepare Output Directory
mkdir(PLOT_DIR_NAME);

log_fid = fopen(PLOT_DIR_NAME + "log.txt","w");


%% Read in data from json-formatted source
covid_file = './covid.json';

TimeSeries = TimeSeriesFromODFile(covid_file,CountriesToTrack);

%% Plot Time Series in Various Ways

conf_case_fig = figure();
set(conf_case_fig,'Name','Confirmed Cases');
legend_for_plot = string([]);
for ii = 1:length(TimeSeries)
    if(MIN_POP_FOR_JOINED_PLOTS < TimeSeries(ii).population)
        legend_for_plot(end+1) = TimeSeries(ii).name;
        semilogy(...
            TimeSeries(ii).dates,...
            TimeSeries(ii).cum_cases...
            );
        hold on;
        text(...
            TimeSeries(ii).dates(end),...
            TimeSeries(ii).cum_cases(end),...
            TimeSeries(ii).name,...
            "Interpreter","none");
    end
end
legend(legend_for_plot,'Location','northwest',"Interpreter","none");
title("Cumulative Cases by Date, Minimum Population " +...
    num2str(MIN_POP_FOR_JOINED_PLOTS),"Interpreter","none");
xlabel("Date");
ylabel("Cumulative Confirmed Cases");
grid on;
saveas(conf_case_fig, PLOT_DIR_NAME + "ConfCases.png");

conf_death_fig = figure();
set(conf_death_fig,'Name','Confirmed Deaths');
legend_for_plot = string([]);
for ii = 1:length(TimeSeries)
    if(MIN_POP_FOR_JOINED_PLOTS < TimeSeries(ii).population)
        legend_for_plot(end+1) = TimeSeries(ii).name;
        semilogy(...
            TimeSeries(ii).dates,...
            TimeSeries(ii).cum_deaths...
            );
        hold on;
        text(...
            TimeSeries(ii).dates(end),...
            TimeSeries(ii).cum_deaths(end),...
            TimeSeries(ii).name,...
            "Interpreter","none");
    end
end
legend(legend_for_plot,'Location','northwest');
title("Cumulative Deaths by Date, Minimum Population " +...
    num2str(MIN_POP_FOR_JOINED_PLOTS),"Interpreter","none");
xlabel("Date");
ylabel("Cumulative Confirmed Deaths");
grid on;
saveas(conf_death_fig, PLOT_DIR_NAME + "ConfDeaths.png");

since_thresh_fig = figure();
set(since_thresh_fig,'Name',"Cases After " + num2str(MIN_INF_THRESH));
legend_for_plot = string([]);
for ii = 1:length(TimeSeries)
    post_thresh_mask = (MIN_INF_THRESH <= TimeSeries(ii).cum_cases);
    if(~any(post_thresh_mask)),continue;end
    if(TimeSeries(ii).population < MIN_POP_FOR_JOINED_PLOTS),continue;end
    legend_for_plot(end+1) = TimeSeries(ii).name;
    thresh_date = TimeSeries(ii).dates(find(post_thresh_mask,1));
    days_since_thresh = days(...
        TimeSeries(ii).dates(post_thresh_mask)-thresh_date...
        );
    semilogy(...
        days_since_thresh,...
        TimeSeries(ii).cum_cases(post_thresh_mask)...
        );
    hold on;
    text(...
        days_since_thresh(end),...
        TimeSeries(ii).cum_cases(end),...
        TimeSeries(ii).name,...
        "Interpreter","none");
end
legend(legend_for_plot,'Location','southeast',"Interpreter","none");
title("Cumulative Cases After Reaching " +...
    num2str(MIN_INF_THRESH) + " Cases, Minimum Population " +...
    num2str(MIN_POP_FOR_JOINED_PLOTS),...
    "Interpreter","none");
xlabel("Days Since " + num2str(MIN_INF_THRESH) + " Cases");
ylabel("Cumulative Confirmed Cases");
grid on;
saveas(since_thresh_fig, PLOT_DIR_NAME + "SinceThresh.png");

for fit_index = 1:length(CountriesToTrack)
    country_to_fit = CountriesToTrack(fit_index);
    cum_cases_here = TimeSeries(fit_index).cum_cases.';
    post_thresh_mask = (MIN_INF_THRESH <= cum_cases_here);
    
    dx = datenum(TimeSeries(fit_index).dates(post_thresh_mask)).';
    dy = TimeSeries(fit_index).cum_deaths(post_thresh_mask).';
    if(isempty(dy)),continue;end
    CF = @(pp,dd) CovidCostFunctional(pp,dx,dy,dd);
    
    p0_assuming_end = [max(dy);0.2;dx(find(max(dy)/2<dy,1))];
    p0_assuming_half = [2*max(dy);0.2;dx(end)];
    p0_choices = [p0_assuming_end,p0_assuming_half];
    
    param_results = 0*p0_choices;
    costs = zeros(1,size(p0_choices,2));
    for p_idx = 1:size(p0_choices,2)
        p0 = p0_choices(:,p_idx);
        [params,nf] = A2CG(CF,p0,0,1e-6,2000);
        param_results(:,p_idx) = params;
        costs(p_idx) = CF(params,0);
    end
    [min_cost,min_idx] = min(costs);
    p0 = p0_choices(:,min_idx);
    params = param_results(:,min_idx);
    
    
    test_alphas = linspace(...
        params(1)/LIN_SEARCH_SCALE,...
        params(1)*LIN_SEARCH_SCALE,...
        LIN_SEARCH_PTS_PER_DIM)';
    test_ks = linspace(...
        params(2)/LIN_SEARCH_SCALE,...
        params(2)*LIN_SEARCH_SCALE,...
        LIN_SEARCH_PTS_PER_DIM)';
    test_hs = linspace(...
        min(dx),...
        min(dx) + LIN_SEARCH_SCALE*(max(dx)-min(dx)),...
        LIN_SEARCH_PTS_PER_DIM)';
    
    [TA,TK,TH] = meshgrid(test_alphas,test_ks,test_hs);
    TC = 0*TA;
    for rr = 1:size(TC,1)
        for cc = 1:size(TC,2)
            for ss = 1:size(TC,3)
                pp = [TA(rr,cc,ss);TK(rr,cc,ss);TH(rr,cc,ss)];
                TC(rr,cc,ss) = CF(pp,0);
            end
        end
    end
    [search_min_cost,l_idx] = min(TC(:));
    [r_idx,c_idx,s_idx] = ind2sub(size(TC),l_idx);    
    
    search_params = [
        TA(r_idx,c_idx,s_idx);
        TK(r_idx,c_idx,s_idx);
        TH(r_idx,c_idx,s_idx)
        ];
    
    [opt_params,~] = A2CG(CF,search_params,0,1e-6,2000);
    
    opt_cost = CF(opt_params,0);
    fprintf("Cost before search: %e\n",min_cost);
    fprintf("Cost after search : %e\n",search_min_cost);
    fprintf("Cost after refine : %e\n\n",opt_cost);
    
    params = opt_params; 
    params_alpha = params(1);
    params_k = params(2);
    params_h = params(3);
    
    Sr = CF(params,0);
    St = sum((dy-mean(dy)).^2);
    r2 = (St-Sr)/St;
    
    
    peak_death_rate = params_k*params_alpha/4;
    fprintf(log_fid,"COVID-19 Predictions for %s\n", country_to_fit);
    fprintf(log_fid,"R-Squared          : %e\n", r2);
    fprintf(log_fid,"Predicted deaths   : %e\n", params(1));
    fprintf(log_fid,"As portion of pop. : %e\n", ...
        params(1)/TimeSeries(fit_index).population);
    fprintf(log_fid,"Day of most deaths : %s\n", datestr(params_h));
    fprintf(log_fid,"Peak death rate    : %e\tdeaths per day\n\n",...
        peak_death_rate);
    
    ff = figure();
    set(ff,'Name',country_to_fit);
    model_dx_min = min(dx);
    p_time = params_h - log(1/FINAL_PLOT_P - 1)/params_k;
    model_dx_max = max(max(dx),p_time);
    model_dx = linspace(model_dx_min,model_dx_max,length(dx));
    dx_as_dates = datetime(dx,'ConvertFrom','datenum');
    model_dx_as_dates = datetime(model_dx,'ConvertFrom','datenum');
    dx_inflection = datetime([params_h,params_h],'ConvertFrom','datenum');
    dy_inflection = [0,params_alpha/2];
    dx_final = datetime([min(model_dx),max(model_dx)],...
        'ConvertFrom','datenum');
    dy_final = [params_alpha,params_alpha];
    plot(dx_as_dates,dy);
    hold on;
    plot(model_dx_as_dates,CovidModel(params,model_dx,0));
    plot(dx_inflection,dy_inflection,'r-.');
    text(dx_inflection(2),dy_inflection(2),...
        "Inflection Point" + newline +...
        datestr(dx_inflection(2)) + newline +...
        num2str(peak_death_rate) + " deaths per day");
    plot(dx_final,dy_final,'r-.');
    text(mean(dx_final),mean(dy_final),...
        newline + "Predicted Count: " +num2str(params_alpha) +...
        " Total Deaths");
    legend("Data","Prediction");
    title("COVID-19 Death Projection for " +country_to_fit,...
        "Interpreter","none");
    xlabel("Date");
    ylabel("Cumulative Deaths");
    grid on;
    saveas(ff,PLOT_DIR_NAME + country_to_fit + ".png");
end

[success,msg,msgid] = rmdir(LATEST_DIR_NAME,'s');
mkdir(LATEST_DIR_NAME);
copyfile(PLOT_DIR_NAME,LATEST_DIR_NAME);

fclose(log_fid);

function M = CovidModel(parameters,dx,degree)
A = parameters(1);
K = parameters(2);
T = parameters(3);
switch(degree)
    case 0
        M = A./(1+exp(-K*(dx-T)));
    case 1
        exp_KT = exp(-K*(dx-T));
        sigm_denom = 1+exp_KT;
        sigm_denom_sq = sigm_denom.^2;
        dMdA = 1./sigm_denom;
        dMdK = A./sigm_denom_sq.*exp_KT.*(dx-T);
        dMdT = -K*A.*exp_KT./sigm_denom_sq;
        M = [dMdA,dMdK,dMdT];
end
end

function C = CovidCostFunctional(parameters,dx,dy,degree)
switch(degree)
    case 0
        % Normal cost functional
        ff = CovidModel(parameters,dx,0);
        rr = dy-ff;
        C = rr'*rr;
    case 1
        % Gradient
        J = CovidModel(parameters,dx,1);
        ff = CovidModel(parameters,dx,0);
        rr = dy-ff;
        C = -2*J'*rr;
end

end