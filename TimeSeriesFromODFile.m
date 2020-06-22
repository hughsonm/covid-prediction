function TimeSeries = TimeSeriesFromODFile(filename,countries_to_track)
covid_text = fileread(filename);
C = jsondecode(covid_text);

for cc = 1:length(countries_to_track)
    ts.name = countries_to_track(cc);
    ts.dates = NaT(0,0);
    ts.new_cases = [];
    ts.new_deaths = [];
    TimeSeries(cc) = ts;
end

current_country = "";
current_country_idx = 0;

for ii = 1:length(C.records)
    iRec = C.records(ii);
    record_country = string(iRec.countriesAndTerritories);
    
    if(record_country ~= current_country)
        current_country = record_country;        
        tracking_this_country = any(contains(countries_to_track,current_country));
        if(tracking_this_country)
            current_country_idx = find(countries_to_track == current_country,1);
            TimeSeries(current_country_idx).population = iRec.popData2019;
        end
    end    
    if(tracking_this_country)
        rec_case_num = (iRec.cases);
        rec_death_num = (iRec.deaths);
        rec_date = datetime(iRec.dateRep,'InputFormat',"dd/MM/yyyy");
        TimeSeries(current_country_idx).dates(end+1) = rec_date;
        TimeSeries(current_country_idx).new_cases(end+1) = rec_case_num;
        TimeSeries(current_country_idx).new_deaths(end+1) = rec_death_num;    
    end    
end

for ii = 1:length(TimeSeries)
    [TimeSeries(ii).dates,sort_idcs] = sort(TimeSeries(ii).dates);
    TimeSeries(ii).new_cases = TimeSeries(ii).new_cases(sort_idcs);
    TimeSeries(ii).new_deaths = TimeSeries(ii).new_deaths(sort_idcs);
    TimeSeries(ii).cum_cases = cumsum(TimeSeries(ii).new_cases);
    TimeSeries(ii).cum_deaths = cumsum(TimeSeries(ii).new_deaths);
    switch(TimeSeries(ii).name)
        case "China"
            % China experienced a fals spike in reported deaths due to a
            % re-classification of some deaths that occurred away from
            % urban hospitals. 
            TS_China = TimeSeries(ii);           
            TS_China.name = "China_Adjusted";
            barnum = 4000;
            dd = TS_China.cum_deaths;
            post_spike_idx = find(barnum<dd,1,'first');
            pre_spike_idx = find(dd<barnum,1,'last');
            scale = dd(post_spike_idx)/dd(pre_spike_idx);
            dd(1:pre_spike_idx) = dd(1:pre_spike_idx)*scale;
            TS_China.cum_deaths = dd;
            TimeSeries(end+1) = TS_China;
    end
    
end
end

