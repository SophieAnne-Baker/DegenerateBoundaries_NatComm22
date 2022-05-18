close all;
clear all;

files = {'3N_curve_data.mat', '5N_curve_data.mat', '7N_curve_data.mat', '9N_curve_data.mat'};
%files = {'curve_data_3N.mat', 'curve_data_5N.mat', 'curve_data_7N.mat', 'curve_data_9N.mat'};

figure
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

hold on
wcses = permute(repmat(reward{1}.wcs, 1, 100, 100), [2, 3, 1]);
scatter(reshape(err{k}, 1, []), reshape(rT{k}, 1, []),15, reshape(wcses, 1, []),'filled')

hold off
clear r r_std reward v v_std dist err_c rT_c
end

ylabel('Mean Decision Time')
xlabel('Mean Error')
title('       SAT curves with {\it N}')
ylim([-2 25])
t = text(0.3, -0.5, '{\it N} =    3,      5,  7, 9');

set(gcf, 'color', 'w', 'units', 'inches', 'position', [0 0 3.75 3.125])

hold on;

ha=gca;

hc = colorbar('Location','eastoutside');
hc.Label.String = '{\it c/W}';
ha.Position = ha.Position - [0 0 .1 .1];
hc.Position = hc.Position + [.001 0 0 0];

if ~exist('figs', 'dir')
       mkdir('figs')
    end

export_fig(append('./figs/', mfilename), '-pdf', '-eps', '-q101');
savefig([pwd '/figs/' mfilename])