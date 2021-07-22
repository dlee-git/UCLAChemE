
function chunk = request_chunk(start_time,end_time,equipment_path,max_count,options,chunk_num)
% URL to get back actual data sent to Academic Hub
% Without maxCount, only 1,000 values per attribute are returned
tic
message = sprintf('Getting chunk %d of 3...', chunk_num);
disp(message);
base_url = 'https://academicpi.azure-api.net/hub/api/';
recorded_url = [ base_url, ...
    'Csv/ElementRecorded',... % controller
    '?path=\\PIAF-ACAD\Classroom Data\Source Data\', equipment_path, ...
    '&startTime=', start_time, '&endTime=', end_time, '&maxCount=', max_count ];
try
chunk = webread(recorded_url, options);  % or interpolated_url 
catch issue
    warning(issue.message)
    return
end
toc
end
