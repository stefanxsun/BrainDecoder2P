%% --- Your preprocessing part (unchanged) ---
mouseTotal = size(data,1);
dayTotal = size(data,2);
cellSel = centerResponsive;
cellSel = visResIdx;

avg_zscore_traces_all = cell(1, dayTotal);

for d = 1:dayTotal
    avg_zscores_day = [];
    for m = 1:mouseTotal
        selected_neurons = cellSel{m};
        selected_neurons = selected_neurons(selected_neurons<=pairedCellNum(m));
        if ~isempty(selected_neurons)
            num_selected_neurons = length(selected_neurons);
            avg_zscores_mouse = NaN(num_selected_neurons, 177);
            for i = 1:num_selected_neurons
                neuron_id = selected_neurons(i);
                zscore_traces = data(m, d).CaCell.Zscore{2, neuron_id};
                if size(zscore_traces,1)>177
                    zscore_traces(178:end,:) = [];
                end
                avg_zscore = mean(zscore_traces, 2, 'omitnan');
                avg_zscores_mouse(i, :) = avg_zscore';
            end
            avg_zscores_day = [avg_zscores_day; avg_zscores_mouse];
        end
    end
    avg_zscore_traces_all{d} = avg_zscores_day;
end

%% --- Visualization ---
figure('Units', 'normalized', 'Position', [0.05 0.1 0.9 0.8]);

% Gather all values across days + diff for consistent color scaling
AllZscoreTraces = [];
for daySel = 1:dayTotal
    AllZscoreTraces = [AllZscoreTraces; avg_zscore_traces_all{daySel}];
end

low = prctile(AllZscoreTraces(:), 1);
high = prctile(AllZscoreTraces(:), 99);

% Layout: upper row = heatmaps, bottom row = average traces
t = tiledlayout(2, dayTotal+1, 'TileSpacing','compact','Padding','compact');

avg_zscores_day1 = avg_zscore_traces_all{1};
avg_zscores_dayN = avg_zscore_traces_all{dayTotal};
avg_zscores_diff = avg_zscores_dayN - avg_zscores_day1;
mean_response = mean(avg_zscores_diff(:, 61:120), 2);
% [~, sort_idx] = sort(mean_response, 'descend');  %%% Use this line to sort for the diff

% --- Heatmaps ---
for daySel = 1:dayTotal
    avg_zscores_day = avg_zscore_traces_all{daySel};
    nexttile(t, daySel); % row 1
    if ~isempty(avg_zscores_day)
        num_cells = size(avg_zscores_day, 1);
        if daySel == 1
            mean_response = mean(avg_zscores_day(:, 61:120), 2);
            [~, sort_idx] = sort(mean_response, 'descend');  %%% Use this line to sort for the first day
        end
        sorted_zscores = avg_zscores_day(sort_idx, :);
        imagesc(sorted_zscores, [0 high]);
        colormap(modified_cmap);
        hold on;
        xline(61, 'w--', 'LineWidth', 1);
        xline(121, 'w--', 'LineWidth', 1);
        hold off;
        title(['Day ' num2str(daySel)]);
        xticks([1, 61, 121, 177]);
        xticklabels({'-2','0','2','4'});
        xlabel('Time (s)');
        ylabel('Cells');
    else
        title(['Day ' num2str(daySel) ' (No Data)']);
    end
end
% Update colorbar (same as before)
cb = colorbar;
cb.Layout.Tile = 'west';
cb.Label.String = 'Average Z-Score';
cb.Label.FontSize = 12;

%% --- Diff heatmap (blue-white-red with gamma stretch) ---
nexttile(t, dayTotal+1);
avg_zscores_day1 = avg_zscore_traces_all{1};
avg_zscores_dayN = avg_zscore_traces_all{dayTotal};
avg_zscores_diff = avg_zscores_dayN - avg_zscores_day1;
sorted_zscores = avg_zscores_diff(sort_idx, :);

imagesc(sorted_zscores);

% --- Symmetric gamma-stretched blue-white-red colormap ---
n = 256;        % total colormap size
gamma = 0.45;    % <1 emphasizes small absolute values

% Symmetric indices from -1 to 1 (0 is center)
idx = linspace(-1, 1, n)';  % -1 -> min, 0 -> zero, 1 -> max

% Apply symmetric gamma: preserve sign
idx_stretch = sign(idx) .* abs(idx).^gamma;

% Interpolate colors
% Linear blue->white->red
blue = [0 0 1];    % min value color
white = [1 1 1];   % zero
red = [1 0 0];     % max value color

cmap_stretched = zeros(n,3);
for i = 1:n
    if idx_stretch(i) < 0
        % negative side: interpolate between blue and white
        mapt = (idx_stretch(i)+1)/1;   % map [-1,0] -> [0,1]
        cmap_stretched(i,:) = blue*(1-mapt) + white*mapt;
    else
        % positive side: interpolate between white and red
        mapt = idx_stretch(i)/1;       % map [0,1] -> [0,1]
        cmap_stretched(i,:) = white*(1-mapt) + red*mapt;
    end
end

% Apply colormap
colormap(gca, cmap_stretched);

% Symmetric color limits
maxVal = max(abs(sorted_zscores(:)));
caxis([-maxVal maxVal]);  % zero is white, symmetric extremes


hold on;
xline(61, 'w--', 'LineWidth', 1);
xline(121, 'w--', 'LineWidth', 1);
hold off;

title('Diff');
xticks([1, 61, 121, 177]);
xticklabels({'-2','0','2','4'});
xlabel('Time (s)');
ylabel('Cells');

% Update colorbar (same as before)
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Diff Day 8 - Day 1';
cb.Label.FontSize = 12;

%% --- Bottom row: Average traces ---
avgTraces = zeros(dayTotal, 177);
for d = 1:dayTotal
    avgTraces(d,:) = mean(avg_zscore_traces_all{d},1,'omitnan');
end
avgDiff = mean(avg_zscores_dayN,1,'omitnan') - mean(avg_zscores_day1,1,'omitnan');

ymin = min([avgTraces(:); avgDiff(:)]);
ymax = max([avgTraces(:); avgDiff(:)]);

for d = 1:dayTotal
    ax = nexttile(t, dayTotal+1+d); % row 2
    plot(avgTraces(d,:), 'LineWidth',1.2); hold on;
    xline(61,'k--'); xline(121,'k--');
    ylim([ymin ymax]);
    title(['Day ' num2str(d)]);
    xticks([1, 61, 121, 177]);
    xticklabels({'-2','0','2','4'});
    xlabel('Time (s)');
    ylabel('Avg resp');
end

% Diff average
nexttile(t, 2*(dayTotal+1));
plot(avgDiff,'k','LineWidth',2); hold on;
xline(61,'k--'); xline(121,'k--');
ylim([ymin ymax]);
title('Diff');
xticks([1, 61, 121, 177]);
xticklabels({'-2','0','2','4'});
xlabel('Time (s)');
ylabel('Avg resp');

sgtitle('Average Z-Score Traces Across Days', 'FontSize', 16);
