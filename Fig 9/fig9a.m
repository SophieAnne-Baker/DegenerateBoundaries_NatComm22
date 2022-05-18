close all;
clear all;

files = {'3N_curve_data.mat', '5N_curve_data.mat', '7N_curve_data.mat', '9N_curve_data.mat'};

C = [100, 143, 255; 220, 38, 127; 255, 176, 0]./255;
for k=1:4
load(files{k})

for wcs = 1:1:length(reward{1}.wcs)
r_temp = r{1}(:,:,wcs);
r_temp(r_temp<-1) = nan;
r_std_temp = r_std{1}(:,:,wcs);
r_std_temp(r_std_temp<-1) = nan;
del = 0.02.*r_std_temp;
r_temp(r_temp < max(max(r_temp)) - del) = nan;
r_opt(:,:,wcs) = r_temp;
end
% get error and RTs for each optimal value
err_c = repmat(squeeze(v{1}(:,:,1)), [1,1,length(reward{1}.wcs)]);
rT_c = repmat(squeeze(v{1}(:,:,2)), [1,1,length(reward{1}.wcs)]);
err_c(isnan(r_opt))=nan;
rT_c(isnan(r_opt))=nan;

err{k} = err_c;
rT{k} = rT_c;

rT_50{k} = reshape(rT{k}(:,:,50), 1, []);
rT_200{k} = reshape(rT{k}(:,:,200), 1, []);
rT_480{k} = reshape(rT{k}(:,:,480), 1, []);

clear r r_std reward v v_std dist err_c rT_c
end

figure; hold on;

logs_k = [log(3+1), log(5+1), log(7+1), log(9+1)]';
rT_50s = [nanmean(rT_50{1}), nanmean(rT_50{2}), nanmean(rT_50{3}), nanmean(rT_50{4})]';
rT_200s = [nanmean(rT_200{1}), nanmean(rT_200{2}), nanmean(rT_200{3}), nanmean(rT_200{4})]';
rT_480s = [nanmean(rT_480{1}), nanmean(rT_480{2}), nanmean(rT_480{3}), nanmean(rT_480{4})]';

for k=1:4
    err_50(k,:)=[min(rT_50{k}), max(rT_50{k})]-rT_50s(k);
    err_200(k,:)=[min(rT_200{k}), max(rT_200{k})]-rT_200s(k);
    err_480(k,:)=[min(rT_480{k}), max(rT_480{k})]-rT_480s(k);
end

p{1} = errorbar(logs_k, rT_50s, err_50(:,1), err_50(:,2), 'color', C(1,:), 'LineWidth', 1);
p{2} = errorbar(logs_k, rT_200s, err_200(:,1), err_200(:,2), 'color', C(2,:), 'LineWidth', 1);
p{3} = errorbar(logs_k, rT_480s, err_480(:,1), err_480(:,2), 'color', C(3,:), 'LineWidth', 1);

for k=1:3
    p{k}.Marker = '.';
    p{k}.MarkerSize = 10;
end

ylabel('Mean Decision Time')
xlabel('log({\it N}+1)')
title('Hick''s Law')
xlim([1.3, 2.4]);

l=legend([p{1}, p{2}, p{3}], {'{\it c/W} = 0.001', '{\it c/W} = 0.04', '{\it c/W} = 0.1'}, 'Location', 'northwest', 'FontSize', 6);

set(gcf, 'color', 'w', 'units', 'inches', 'position', [0 0 2.75 2.5].*1.25)

if ~exist('figs', 'dir')
       mkdir('figs')
    end

export_fig(append('./figs/', mfilename), '-pdf', '-eps', '-q101');
savefig([pwd '/figs/' mfilename])