function [] = MakeLatestMarkdown(location)
%MAKELATESTMARKDOWN Makes a new README.md file display all the svg in a loc

fstr = "";

ds = dir(location);

url_prefix = ...
    "https://raw.githubusercontent.com/"+...
    "hughsonm/covid-prediction/master/Latest/";

for ii = 1:length(ds)
    if ~(ds(ii).name == "." || ds(ii).name == "..")
        [~,country_name,ext] = fileparts(ds(ii).name);
        if ext == ".svg"
            istr = "!["+country_name+"](" + url_prefix + ds(ii).name + ")";
            fstr = fstr + istr + newline;
        end
    end
end

fid = fopen(fullfile(location,"README.md"),"w");
fprintf(fid,"%s",fstr);
fclose(fid);


end

