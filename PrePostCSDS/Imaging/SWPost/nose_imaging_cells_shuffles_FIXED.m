function nose_imaging_cells_shuffles_FIXED(alldeltaF_social,time_stamps,mouse,shuffle_time,prct)
% Find average F signal during attack events across all days 
disp('Running shuffling code...')

%% Z-score
F_mean = nanmean(alldeltaF_social.C_raw_sync,2);
F_std = nanstd(alldeltaF_social.C_raw_sync,0,2);

for numcell = 1:size(alldeltaF_social.C_raw_sync,1)
    alldeltaF_social.C_raw_zscore(numcell,:) = (alldeltaF_social.C_raw_sync(numcell,:)-F_mean(numcell,:))/F_std(numcell,:);
end

%% Time and FR
time = 0:1/100:size(alldeltaF_social.C_raw_zscore,2)/100; % real frame rate is 100fps
time = time(1:end-1);

fps = 30;  % frame rate of the video is 100fps for Motif
shuffle_frames = shuffle_time*fps; % convert shuffle time to frames

% fix any timestamp issues
if length(time_stamps)==1 
    time_stamps(2) = size(alldeltaF_social.C_raw_zscore,2);
elseif rem(length(time_stamps), 2) ~= 0
    time_stamps(end+1) = size(alldeltaF_social.C_raw_zscore,2);
end

%% Find number of steps for 1 sec
unit = fps; %time(2) - time(1);
steps_05 = round(fps)/2;
steps_1 = round(fps);%floor(1/unit);
steps_2 = steps_1*2;%floor(2/unit);
steps_3 = steps_1*3;%floor(3/unit);
steps_5 = steps_1*5;%floor(5/unit);
steps_10 = steps_1*10;%floor(5/unit);
steps_30 = steps_1*30;%floor(30/unit);
steps_02 = steps_1/5;

if ~isempty(time_stamps) % avoid errors due to no attack
    %% Find each behavioral epoch 
    startstop = [round(time_stamps(1:2:end)*fps),round(time_stamps(2:2:end)*fps)];
    
    %% Merge attack sessions that are within 5s of each other 
    k = 1;
    startstop_merge(1,:) = startstop(1,:);
    for j = 2:size(startstop,1)
        if startstop(j,1) - startstop_merge(k,2) < steps_5 % if the next attack start earlier than 5s before the end of the previous attack, merge them into one attack events
            startstop_merge(k,2) = startstop(j,2);
        else
            k = k + 1;
            startstop_merge(k,:) = startstop(j,:);
        end
    end

    %% Find indices on F vector 
    %startstop_F(:,1) = startstop(:,1);
    %startstop_F(:,2) = startstop(:,2);
    startstop_F(:,1) = startstop_merge(:,1);
    startstop_F(:,2) = startstop_merge(:,2);
    
    if startstop_F(end,2) > size(alldeltaF_social.C_raw_zscore,2)
        startstop_F(end,2) = size(alldeltaF_social.C_raw_zscore,2);
    end
    
    % event start
    %startstop_F_start(:,1) = startstop(:,1);
    startstop_F_start(:,1) = startstop_merge(:,1)-steps_2; % steps_5; % take 2 sec beofre start of all attacks 
    startstop_F_start(:,2) = startstop_merge(:,1)+steps_5; % steps_5; % take 5 sec after start of all attacks 

    % event end
    %startstop_F_end(:,1)= startstop(:,2);
    startstop_F_end(:,1)= startstop_merge(:,2)-steps_2; % steps_5; % take 2 sec before end of all attacks 
    startstop_F_end(:,2) = startstop_merge(:,2)+steps_5; % steps_5; % take 5 sec after end of all attacks 
    
    % avoid errors due to behavior starting before session %< 1 sec after start of session
    if startstop_F_start(1,1) <= 0 %steps_1
        startstop_F_start(1,:) = [];
    end

    if startstop_F_end(1,1) <= 0 %steps_2
        startstop_F_end(1,:) = [];
    end
    
    % avoid errors due to behavior starting < 1 sec from end of session
    if ~isempty(startstop_F_start)
        if startstop_F_start(end,2) > size(alldeltaF_social.C_raw_zscore,2)
            startstop_F_start(end,:) = [];  
        end
    end

    if ~isempty(startstop_F_end)
        if startstop_F_end(end,2) > size(alldeltaF_social.C_raw_zscore,2)
            startstop_F_end(end,:) = [];
        end
    end

   %% Create a matrix of responses to event start 
    clear start_e;
    clear stop_e;
    start_e(:,1) = startstop_F_start(:,1);  % 2 sec prior to event start
    stop_e(:,1) = startstop_F_start(:,2); % 5 sec post event start
    
    % error checking: make sure event doesn't start before trial starts
    err1 = find(start_e<=0);

    if ~isempty(err1)
        start_e(err1,:) = [];
        stop_e(err1,:) = [];
        startstop_F_start(err1,:) = [];
    end
    
    err2 = find(stop_e>length(time));

    if ~isempty(err2)
        start_e(err2,:) = [];
        stop_e(err2,:) = [];
    end
      
    if sum(~isnan(startstop_F_start)) > 1
        % cell bodies
        for j = 1:size(start_e,1)
            for cells = 1:size(alldeltaF_social.C_raw_zscore,1)
                avg_start(cells,j,:) = alldeltaF_social.C_raw_zscore(cells,start_e(j):stop_e(j));
            end
        end
    else
        avg_start = nan(size(start_e,1),steps_5*2+1);%steps_1*2+1);
    end

    %% Find start of attack Averages
    % cell bodies
    avg_start_all = squeeze(nanmean(avg_start,2));
    std_start_all = squeeze(nanstd(avg_start,0,2)/sqrt(size(avg_start,2)));

    %% create 1000 shuffles for each cell
    C_raw_shuffles = circshift(alldeltaF_social.C_raw_zscore,shuffle_time*fps,2);
    for s = 2:1000
        C_raw_shuffles(:,:,s) = circshift(C_raw_shuffles(:,:,s-1),shuffle_time*fps,2);
    end

    %% Plot some of the shuffles
    % % m569
    % figure('DefaultAxesFontSize',12)
    % pos8 = [0.2 0.54 0.3 0.1];
    % subplot('Position',pos8)
    % numcell = 5; 
    % time_len = size(C_raw_shuffles,2);
    % time_plot = 0:fps:time_len/fps;
    % imagesc(squeeze(C_raw_shuffles(numcell,:,1:10))')
    % colorbar
    % clim([-2 2])
    % a=colorbar;
    % xticks(0:fps*30:time_len)
    % xticklabels(time_plot)
    % ylabel(a,'\DeltaF/F (Z)','FontSize',12,'Rotation',90);
    % xlabel('Time (s)')
    % ylabel('Shuffles')
    % title('Example cell (shifted 5s)')
    % box off
    
    %% Get avg trace for all shuffles across all timestamps
     if sum(~isnan(startstop_F_start)) > 1
        % cell bodies
        for s = 1:size(C_raw_shuffles,3)
            for j = 1:size(start_e,1)
                for cells = 1:size(alldeltaF_social.C_raw_zscore,1)
                    avg_start_shuffle{s,1}(cells,j,:) = C_raw_shuffles(cells,start_e(j):stop_e(j),s);
                end
            end
        end
    else
        avg_start_shuffle = nan(size(start_e,1),steps_5*2+1);%steps_1*2+1);
     end

    %% Find shuffled start of attack Averages
    % cell bodies
    for s = 1:size(C_raw_shuffles,3)
        avg_start_shuffle_all(:,:,s) = squeeze(nanmean(avg_start_shuffle{s},2));
    end

    %% Plot timelocked traces
    % numcell = 5; %m569
    % time_plot = -2:2:5;
    % figure('DefaultAxesFontSize',12)
    % pos8 = [0.2 0.54 0.1 0.1];
    % subplot('Position',pos8)
    % for s = 1:100%size(C_raw_shuffles,3)
    %     plot(squeeze(avg_start_shuffle_all(numcell,:,s)),'Color',[.7 .7 .7],'LineWidth',0.25)
    %     hold on
    % end
    % plot(avg_start_all(numcell,:),'r','LineWidth',1)
    % plot([60,60],[-2.5 4.5],'k--','linewidth',1)
    % xlim([0 210])
    % ylim([-2 2.5])
    % xticks(0:60:210)
    % %yticks(-0.5:.5:1.5)
    % xticklabels(time_plot)
    % xlabel('Time from social zone entry (s)')
    % ylabel('Mean \DeltaF/F (Z)')
    % title('Example Cell vs 100 shuffles')
    % box off

    %% Find mean value after attack onset for each cell and all shuffles (start)
    start_p = steps_2+steps_05;%steps_1;  % start 1s after attack start
    stop_p = start_p+steps_2;%steps_1; % stop 2s after attack start

    % subtract avg of pre-attack before finding peak
    % avg_start_all_sub = avg_start_all - nanmean(avg_start_all(:,1:steps_5),2);

    % find location of the real peak for each cell 1-2s after attack onset
    % start
    pk_real = nanmean(avg_start_all(:,start_p:stop_p),2);

    % find the peak for each shuffle of the data
    pk_shuffles = squeeze(nanmean(avg_start_shuffle_all(:,start_p:stop_p,:),2));

    %% Find the 95 percentile of the real peak and shuffles vectors
    sig_val_start = prctile([pk_real,pk_shuffles],prct,2);
    sig_up_start = pk_real>sig_val_start;
    for_up_plot = sig_val_start(numcell);

    sig_val_start = prctile([pk_real,pk_shuffles],100-prct,2);
    sig_down_start = pk_real<sig_val_start;
    for_down_plot = sig_val_start(numcell);

    %% Plot a histogram of shuffled peak averages vs real peak average
    % numcell = 5; %4;
    % figure('DefaultAxesFontSize',12)
    % pos8 = [0.2 0.54 0.1 0.1];
    % subplot('Position',pos8)
    % y_val = 150;
    % histogram(pk_shuffles(numcell,:),'FaceColor',[.7 .7 .7])
    % hold on
    % plot([pk_real(numcell),pk_real(numcell)],[0,y_val],'r','LineWidth',1)
    % plot([for_up_plot,for_up_plot],[0,y_val])
    % plot([for_down_plot,for_down_plot],[0,y_val])
    % xlabel('Mean \DeltaF/F (Z)')
    % title('Example cell vs null distribution')
    % box off

    %% Plot 
    % numcell = 3;
    % histogram(pk_shuffles(numcell,:))
    % hold on
    % plot([pk_real(numcell),pk_real(numcell)],[0,160],'k','LineWidth',1)
else
    sig_up_start = nan(size(alldeltaF_social.C_raw_zscore,1),1);
    sig_down_start = nan(size(alldeltaF_social.C_raw_zscore,1),1);
    avg_start = nan(size(alldeltaF_social.C_raw_zscore,1),1,steps_5*2+1);
    avg_start_all = squeeze(nan(size(avg_start)));
    pk_real = nan(size(alldeltaF_social.C_raw_zscore,1),1);
end
%% save
save(strcat(mouse,'_sig_nose_cells_FIXED.mat'),'sig_up_start','sig_down_start','avg_start','avg_start_all',"alldeltaF_social",'pk_real')