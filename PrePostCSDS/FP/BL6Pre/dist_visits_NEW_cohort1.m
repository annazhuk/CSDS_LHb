function dist_visits_NEW_cohort1(dist_center_filtered,F_sync,mouse,strain)

%% time and fps
fps = 29.97; 
steps_1 = fps;
steps_2 = steps_1*2;
steps_5 = steps_1*5;
time = 0:1/fps:300;

%% Normalize and z-score
% Find the mean for a rolling window 
F_mean = movmean(F_sync,1000,'omitnan');

F_final = nan(1,length(F_sync));

for i = 1:length(F_final)
    F_final(i) = (F_sync(i) - F_mean(i))/F_mean(i);
end

% z-score
F_mean = nanmean(F_final);
F_std = nanstd(F_final);

z_score = (F_final - F_mean)/F_std;

%% Fdata and time may be off by a few frames. Make sure they're the same length
if length(z_score)<length(time)
    z_score(length(z_score)+1:length(time)) = nan;
elseif length(z_score)>length(time)
    z_score(length(time)+1:length(z_score)) = [];
end

%% Find start and stop frames when mouse gets with a certain distance to cups
close_frames = find(dist_center_filtered<=9);

if ~isempty(close_frames)
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
    startstop_merge(1,:) = startstop_crop(1,:);
    for j = 2:size(startstop_crop,1)
        if startstop_crop(j,1) - startstop_merge(k,2) < steps_1 % if the next entrance into the close distance start earlier than 1s before the end of the previous entrance, merge them into one
            startstop_merge(k,2) = startstop_crop(j,2);
        else
            k = k + 1;
            startstop_merge(k,:) = startstop_crop(j,:);
        end
    end

    %% Save visits for plot
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
        elseif startstop_plot(i,2) > length(z_score)
            errs(k,1) = i;
            k = k + 1;
        end
    end

    startstop_plot(errs,:) = [];

    %% Save visits into a matrix
    visits_plot = nan(size(startstop_plot,1),startstop_plot(1,2) - startstop_plot(1,1)+1);

    for i = 1:size(startstop_plot,1)
        visits_plot(i,:) = z_score(startstop_plot(i,1):startstop_plot(i,2));
    end

    visits_plot_avg = nanmean(visits_plot,1);
else
    visits_plot_avg = nan(1,round(steps_2+steps_5)+1);
end

%% save
save(strcat(mouse,'_',strain,'_visitsplot_NEW.mat'),'visits_plot_avg','visits_plot')