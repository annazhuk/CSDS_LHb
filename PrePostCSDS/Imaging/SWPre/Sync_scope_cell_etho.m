function Sync_scope_cell_etho(filename,alldeltaF_social,mousename,cond,scope_FR,videoname)
% sync wholefield signal to camera 
%% Import data and get Frame Rate
% Uncropped F data

%TracerDAQ data 
[Sample,camera,scope] = importfile_ImagingDAQ(filename);  

% Import video
video = VideoReader(videoname); 
camera_FR = video.FrameRate; % camera frame rate (can get it from the video file)
%camera_FR = 25; % Motif FR is 25fps

%% Get start/stop times
% Obtain Tracer daq acquisition rate from user
prompt={'Enter Tracer Daq acquisition rate'};
name2 = 'Tracer Daq acquisition rate';
defaultans2 = {'1000'};
options.Interpreter = 'tex';
answer = inputdlg(prompt,name2,[1 40],defaultans2,options);
Acq_rate=str2num(answer{1,1});
%TTL_unit = 1000; %etho_length/video_length;  % 1 second per each TTL unit

t=Sample/Acq_rate; %(time in seconds)
plot(t,camera); 

% Obtain Threshold to detect peaks
prompt={'Enter a threshold'};
name = 'Threshold_camera';
defaultans = {'1'}; %4.78
options.Interpreter = 'tex';
answer = inputdlg(prompt,name,[1 40],defaultans,options);
threshold=str2num(answer{1,1});

close all;

[testx,testy]=find(camera(:,1)>threshold,1,'first');
camera_start_index=(testx-1);
camera_start_time= (testx-1)/Acq_rate;
[testx2,testy2]=find(camera(:,1)>threshold,1,'last');
camera_stop_time= (testx2+1)/Acq_rate;
camera_stop_index=(testx2+1);

plot(t,camera); 
hold on
title('Make sure that the start and stop pulses are being detected correctly');
plot(camera_start_time,camera(round(camera_start_time*Acq_rate),1),'ko')
plot(camera_stop_time,camera(round(camera_stop_time*Acq_rate),1),'ko')

video_length = camera_stop_time - camera_start_time + 0.1;

plot(t,scope)

[testx,testy]=find(scope(:,1)>threshold,1,'first');
scope_start_index=(testx-1);
scope_start_time= (testx-1)/Acq_rate;
[testx2,testy2]=find(scope(:,1)>threshold,1,'last');
scope_stop_time= (testx2+1)/Acq_rate;
scope_stop_index=(testx2+1);

plot(scope_start_time,scope(round(scope_start_time*Acq_rate),1),'ro')
plot(scope_stop_time,scope(round(scope_stop_time*Acq_rate),1),'ro')

%% Error checking: make sure that imaging scope is Ch2 and video camera is Ch3
lag = camera_start_index - scope_start_index;

if lag < 0
    error('Video camera starting before imaging scope. Check channels in tracer DAQ .csv file!')
end

%% Find each camera and scope pulse
etho_delay_start = 0.005; %check ethovision trial control settings for this setting
etho_delay_end = 0.005; %check ethovision trial control settings for this setting
camera_block = Acq_rate*(etho_delay_start+etho_delay_end/2);

scope_delay_start = 0.02;
scope_delay_end = 0.02;
scope_block = Acq_rate*(scope_delay_start+scope_delay_end/2);

[camera_peaks,camera_loc] = findpeaks(camera,'MinPeakHeight',threshold,'MinPeakDistance',(camera_block));
%plot(camera_loc/Acq_rate,camera_peaks,'bo')

[scope_peaks,scope_loc] = findpeaks(scope,'MinPeakHeight',threshold,'MinPeakDistance',(scope_block));
%plot(camera_loc/Acq_rate,camera_peaks,'bo')

%% Error checking 
% make sure number of scope frames is the same as the F vector
if length(scope_loc) ~= length(alldeltaF_social.C_raw(1,:))
    warning('number of scope frames and length of F vector is not the same!')
    diff = length(scope_loc) - length(alldeltaF_social.C_raw(1,:))
end

%% Crop scope to camera
% start frame
temp = scope_loc - camera_start_index;
scope_start_crop = find(temp==min(temp(temp>0)));


% stop frame
temp = camera_stop_index - scope_loc;
scope_stop_crop = find(temp==min(temp(temp>0)));

%% Find the very last pulse
% Find where the video should have started based on duration TTL_unit
delay = 5; % added a 5s delay to solve the partial frame issue
vid_start = camera_stop_index - ((video.Duration-delay)*Acq_rate);
frames = (video.Duration-delay)*camera_FR; % total number of frames

% Find the location of each frame using ending pulse, number of frames, and
% vid_start
y = double(uint64(linspace(camera_start_index,camera_stop_index,frames)))';

%% Error checking
% Check to see how far off video time is from number of pulses found 
figure;
plot(camera)
hold on
plot(camera_loc,camera_peaks,'bo');
plot(y,camera(y),'ro');

if (video.Duration - delay) == video_length
else 
    warning('Incorrect number of pulses based on video length!')
    missing_time = (video.Duration-delay)-video_length
end

%% Interpolate F data to fit with the camera frames
x = scope_loc(scope_start_crop:scope_stop_crop);
xq = y;
v_raw = double(alldeltaF_social.C_raw(:,scope_start_crop:scope_stop_crop));
v_C = double(alldeltaF_social.C(:,scope_start_crop:scope_stop_crop));
v_S = double(alldeltaF_social.S(:,scope_start_crop:scope_stop_crop));

for cell = 1:size(alldeltaF_social.C_raw,1)
    alldeltaF_social.C_raw_sync(cell,:) = interp1(x,v_raw(cell,:),xq,'spline');
    alldeltaF_social.C_sync(cell,:) = interp1(x,v_C(cell,:),xq,'spline');
    alldeltaF_social.S_sync(cell,:) = interp1(x,v_S(cell,:),xq,'spline');
end

%% Name the file
lastpart = 'celldata_';
name = strcat(lastpart,mousename,'_',cond);

%% Save
save (name,'alldeltaF_social','mousename') %,'video');
end
