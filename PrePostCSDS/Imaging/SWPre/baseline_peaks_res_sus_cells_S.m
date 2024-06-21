function baseline_peaks_res_sus_cells_S(allmice_im,SW_CSDS_Post,mouse_CSDS,SW_CTRL_Post,CSDS_ID,CSDS_SI)

%% Find where to split
avg = nanmean(SW_CTRL_Post);
std_val = nanstd(SW_CTRL_Post);
med = avg - std_val; %med = floor(avg - std_val);


%med = median(SW_CSDS_Post);

SI_CSDS = SW_CSDS_Post;

resID = mouse_CSDS(find(SI_CSDS>=med));
susID = mouse_CSDS(find(SI_CSDS<med));

res_im = struct;
sus_im = struct;

%% Split mice into res and sus
% cell bodies
% res
r = 1;
for i = 1:size(allmice_im.mice,1)
    if ismember(allmice_im.mice(i),resID)
        res_im.mice(r,1) = allmice_im.mice(i);
        res_im.pk_totals{r,1} = allmice_im.pk_totals{i,1};
        r = r + 1;
    end
end

% sus
s = 1;
for i = 1:size(allmice_im.mice,1)
    if ismember(allmice_im.mice(i),susID)
        sus_im.mice(s,1) = allmice_im.mice(i);
        sus_im.pk_totals{s,1} = allmice_im.pk_totals{i,1};
        s = s + 1;
    end
end

%% Sort
% Crop out mice npt analyzed (bad signal, gfp, etc)
% load SI data and make sure the IDs and SIs are in the correct order!

CSDS_SI = CSDS_SI/300*100;

% cell bodies
cropped_ID = CSDS_ID(find(ismember(CSDS_ID,allmice_im.mice)));
cropped_SI = CSDS_SI(find(ismember(CSDS_ID,allmice_im.mice)));
% make sure croppedSIs are in the same order as struct with Favg data
[~,idx] = sort(allmice_im.mice);

% Favg data and IDs should be in order, but double check this 
if sum(idx == (1:length(idx))') == length(idx)
    disp('sorted order is correct')
else
    warning('Error in the mouse ID ordering!')
end

[cropped_ID_sorted_im,idx_sorted] = sort(cropped_ID);
cropped_SI_sorted_im = cropped_SI(idx_sorted);

%% save
save('baseline_peaks_res_sus_cells_S_NEW.mat','sus_im','res_im','allmice_im','cropped_ID_sorted_im','cropped_SI_sorted_im')