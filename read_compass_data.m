function [comp] = read_compass_data(fname)
% function to read compass data file

% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 14);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ";";

% Specify column names and types
opts.VariableNames = ["date", "timelocaltime", "altitudem", "longitude", "latitude", "tiltdeg", "rolldeg", "compassdeg", "gyroxdegs", "gyroydegsec", "gyrozdegs", "pressurehPa", "speedms", "accelerationms"];
opts.VariableTypes = ["datetime", "datetime", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "date", "InputFormat", "yyyy-MM-dd");
opts = setvaropts(opts, "timelocaltime", "InputFormat", "HH:mm:ss.SSS");
opts = setvaropts(opts, "pressurehPa", "TrimNonNumeric", true);
opts = setvaropts(opts, "pressurehPa", "ThousandsSeparator", ",");

% Import the data
GPS = readtable(fname, opts);

for i = 1:size(GPS,1)
    comp.time(i) = datenum(join([datestr(GPS{i,1}),' ',datestr(GPS{i,2},'HH:MM:SS.FFF')]));
end
compassdata = table2array(GPS(:,3:end));
comp.lat = compassdata(:,3);
comp.lon = compassdata(:,2);
comp.tilt = compassdata(:,4);
comp.roll = compassdata(:,5);
comp.compass = compassdata(:,6);
