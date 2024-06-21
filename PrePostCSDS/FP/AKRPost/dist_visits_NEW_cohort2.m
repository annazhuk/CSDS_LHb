function dist_visits_NEW_cohort2(dist_center_filtered,F_sync_cb,mouse,strain)

%% time and fps
fps = 29.97; 
steps_1 = fps;
steps_2 = steps_1*2;
steps_5 = steps_1*5;
time = 0:1/fps:300;

%% Normalize and z-score
% Find the mean for a rolling window 
F_mean_cb = movmean(F_sync_cb,1000,'omitnan');

F_final_cb = nan(1,length(F_sync_cb));

for i = 1:length(F_final_cb)
    F_final_cb(i) = (F_sync_cb(i) - F_mean_cb(i))/F_mean_cb(i);
end

% z-score
F_mean_cb = nanmean(F_final_cb);
F_std_cb = nanstd(F_final_cb);

z_score_cb = (F_final_cb - F_mean_cb)/F_std_cb;

%% Fdata and time may be off by a few frames. Make sure they're the same length
if length(z_score_cb)<length(time)
    z_score_cb(length(z_score_cb)+1:length(time)) = nan;
elseif length(z_score_cb)>length(time)
    z_score_cb(length(time)+1:length(z_score_cb)) = [];
end

%% Find start and stop frames when mouse gets with a certain distance to cups
close_frames = find(dist_center_filtered<=9);

idx_stop = close_frames(find(diff(close_frames)>1));
idx_start = close_frames(find(diff(close_frames)>1)+1);
startstop(:,1) = [close_frames(1);idx_start];
startstop(:,2) = [idx_stop;close_frames(end)];

%% Crop out short visits 
visit_lengths = startstop(:,2) - startstop(:,1);
startstop_crop = startstop;
startstop_crop(visit_lengths<steps_2,:) = [];

%% Combine visits that are within a sort time frame of one another 
k = 1;
if ~isempty(startstop_crop)
    startstop_merge(1,:) = startstop_crop(1,:);
    for j = 2:size(startstop_crop,1)
        if startstop_crop(j,1) - startstop_merge(k,2) < steps_1 % if the next entrance into the close distance start earlier than 1s before the end of the previous entrance, merge them into one
            startstop_merge(k,2) = startstop_crop(j,2);
        else
            k = k + 1;
            startstop_merge(k,:) = startstop_crop(j,:);
        end
    end
else
    startstop_merge = startstop_crop;
end

%% Save visits for plot
if ~isempty(startstop_crop)
    startstop_plot(:,1) = startstop_crop(:,1) - round(steps_2);
    startstop_plot(:,2) = startstop_crop(:,1) + round(steps_5);
    
    % avoid errors by getting rid of visits that start right at the beginning of
    % the trial and end right at the end
    errs = [];
    k = 1;
    for i = 1:size(startstop_plot,1)
        if startstop_plot(i,1) <= 0 
            errs(k,1) = i;
            k = k + 1;
        elseif startstop_plot(i,2) > length(z_score_cb)
            errs(k,1) = i;
            k = k + 1;
        end
    end
    startstop_plot(errs,:) = [];
end

%% Save visits into a matrix
if ~isempty(startstop_crop)
    visits_plot = nan(size(startstop_plot,1),startstop_plot(1,2) - startstop_plot(1,1)+1);
    
    for i = 1:size(startstop_plot,1)
        visits_plot(i,:) = z_score_cb(startstop_plot(i,1):startstop_plot(i,2));
    end
    visits_plot_avg = nanmean(visits_plot,1);
else
    visits_plot = [];
    visits_plot_avg = nan(1,round(steps_2)+round(steps_5)+1);
end

%% save
save(strcat(mouse,'_',strain,'_visitsplot_NEW.mat'),'visits_plot_avg','visits_plot')
