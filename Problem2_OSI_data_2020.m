clc,clear
%% Use this script to access lab data from the OSIsoft Academic Hub
% Formatted time-series data will be returned to the table called 'data'
%---------------------------------------------------

%% Parameters to modify for data access 
% For 107 project, please do not change time range and maxcount in this section
% After run the script, please use
% variable named 'data' to calculate mean and std for both open-loop and
% closed-loop data
% Please use the latest Matlab version to run this script. Some older
% versions (e.g., 2015b) may not work.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time can be specified in 2 ways: absolute or relative 
% 1a. Absolute time (50,000 first values)

% 1b. 50,000 last values (uncomment to use)
%StartTime = '2016-01-23T03:59:24Z'; % GMT 
%EndTime = '2016-02-15T00:00:00Z';

% StartTime = '1-Jan-2018 00:00:00+GMT';
% EndTime = '1-Jan-2018 00:20:00+GMT';
% 2. Relative (effective, comment if you use absolute)
%     Suffixex are 'h', 'm' and 's' for hour, minute and second respectively
% StartTime = '*-1h'; % meaning is 'an hour ago' since
% EndTime   = '*';    %  '*' stands for 'now'

%Note on Timestamp format: OSIsoft Academic Hub uses UTC, ISO8601 format 
%For more information on ISO8601: https://academic.osisoft.com/time-help

Interval = '9s'; % '1s' for 1-second interval; '1m' for 1-minute interval
% 25ms for a 25 millisecond interval. *** Do not use interval < 25ms ***

% this is the path to where your sensor data is stored on the Academic Hub
EquipmentPath= 'UCLA\Open_Loop';  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change EquipmentPath to one of the following names, to download three
% sets of data
% 'UCLA\Open_Loop', 'UCLA\Close_Loop_1', 'UCLA\Close_Loop_2'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% those are McMaster data access credentials 
username = 'uclareader0'; 
password = 'qTM8yabzs4$Wtx';

% Controls the maximum value per attribute to return. Failure will occur 
%  if using a value > 50,000 for MaxCount
MaxCount = '20000';

%% Data access portion - students should not need to edit this

options = weboptions('Username', username,...
    'Password', password, 'ContentType', 'auto',...
    'Timeout', 30, 'CertificateFilename', '',...
    'HeaderFields', ["Accept-Encoding" ""]);

StartTime = '2020-01-01T00:00:00'; % Pacific Time Zone
EndTime   = '2020-01-14T07:59:24';
data1 = request_chunk(StartTime,EndTime,EquipmentPath,MaxCount,options,1);
StartTime = '2020-01-14T08:00:00'; % Pacific Time Zone
EndTime   = '2020-01-28T15:59:24';
data2 = request_chunk(StartTime,EndTime,EquipmentPath,MaxCount,options,2);
StartTime = '2020-01-28T16:00:00'; % Pacific Time Zone
EndTime   = '2020-02-11T20:00:00';
data3 = request_chunk(StartTime,EndTime,EquipmentPath,MaxCount,options,3);
%%
tic
disp('Rebuilding CSV...')
data = [data1; data2; data3];
clear data1 data2 data3;
a = data.Attribute;
index =zeros(180000,1);
for k = 1:180000
index(k) = (a{k,1}~= "Heating" );
end
data = data(logical(index),:);
attributes = unique(data.Attribute);
keep_idx = arrayfun(@(x) strcmp(x,attributes(1)), data.Attribute);
timestamps = data(keep_idx, :).Timestamp;
values = zeros(numel(timestamps), numel(attributes));
for i = 1:numel(attributes)
    keep_idx = arrayfun(@(x) strcmp(x,attributes(i)), data.Attribute);
    values(:, i) = data(keep_idx, :).Value;
end
%%
headers = ['Timestamp'; attributes];
final_csv = table(timestamps);
for i = 1:numel(attributes)
    final_csv = addvars(final_csv, values(:, i));
end
final_csv.Properties.VariableNames = arrayfun(@(x) matlab.lang.makeValidName(x), headers);
% final_csv 
toc
