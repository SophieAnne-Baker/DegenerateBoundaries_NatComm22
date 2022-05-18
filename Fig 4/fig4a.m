clear all;
close all;

load('3N_curve_data_dense.mat')
%seprarate optima for each wcs into matrix r_opt
for wcs = 1:1:length(reward.wcs)
    r_temp = r(:,:,wcs);
    r_temp(r_temp<-1) = nan;
    r_std_temp = r_std(:,:,wcs);
    r_std_temp(r_std_temp<-1) = nan;
    del = 0.02.*r_std_temp;
    r_temp(r_temp < max(max(r_temp)) - del) = nan;
    r_opt(:,:,wcs) = r_temp;
end
% get error and RTs for each optimal value
err_c = repmat(squeeze(v(:,:,1)), [1,1,length(reward.wcs)]);
rT_c = repmat(squeeze(v(:,:,2)), [1,1,length(reward.wcs)]);
err_c(isnan(r_opt))=nan;
rT_c(isnan(r_opt))=nan;
clear r_temp
clear r_opt
clear r_std_temp

load('3N_oscil_data.mat')
%seprarate optima for each wcs into matrix r_opt
for wcs = 1:1:length(reward.wcs)
    r_temp = r(:,:,:,wcs);
    r_temp(r_temp<-1) = nan;
    r_std_temp = r_std(:,:,:,wcs);
    r_std_temp(r_std_temp<-1) = nan;
    del = 0.02.*r_std_temp;
    r_temp(r_temp < max(max(r_temp)) - del) = nan;
    r_opt(:,:,:,wcs) = r_temp;
end
% get error and RTs for each optimal value
err_o = repmat(squeeze(v(:,:,:,1)), [1,1,1,length(reward.wcs)]);
rT_o = repmat(squeeze(v(:,:,:,2)), [1,1,1,length(reward.wcs)]);
err_o(isnan(r_opt))=nan;
rT_o(isnan(r_opt))=nan;
clear r_temp
clear r_opt
clear r_std_temp

load('3N_power_data.mat')
%seprarate optima for each wcs into matrix r_opt
for wcs = 1:1:length(reward.wcs)
    r_temp = r(:,:,:,wcs);
    r_temp(r_temp<-1) = nan;
    r_std_temp = r_std(:,:,:,wcs);
    r_std_temp(r_std_temp<-1) = nan;
    del = 0.02.*r_std_temp;
    r_temp(r_temp < max(max(r_temp)) - del) = nan;
    r_opt(:,:,:,wcs) = r_temp;
end
% get error and RTs for each optimal value
err_p = repmat(squeeze(v(:,:,:,1)), [1,1,1,length(reward.wcs)]);
rT_p = repmat(squeeze(v(:,:,:,2)), [1,1,1,length(reward.wcs)]);
err_p(isnan(r_opt))=nan;
rT_p(isnan(r_opt))=nan;
clear r_temp
clear r_opt
clear r_std_temp


for i = 1:length(reward.wcs)
    c_err_av(i) = nanmean(nanmean( err_c(:,:,i) ));
    o_err_av(i) = nanmean(nanmean(nanmean( err_o(:,:,:,i) )));
    p_err_av(i) = nanmean(nanmean(nanmean( err_p(:,:,:,i) )));
    err_av(i) = nanmean([c_err_av(i) o_err_av(i) p_err_av(i)]);
    
    c_rT_av(i) = nanmean(nanmean( rT_c(:,:,i) ));
    o_rT_av(i) = nanmean(nanmean(nanmean( rT_o(:,:,:,i) )));
    p_rT_av(i) = nanmean(nanmean(nanmean( rT_p(:,:,:,i) )));
    rT_av(i) = nanmean([c_rT_av(i) o_rT_av(i) p_rT_av(i)]);
end

figure
hold on
wcses = permute(repmat(reward.wcs, 1, 100, 100), [2, 3, 1]);
scatter(reshape(err_c, 1, []), reshape(rT_c, 1, []),15, reshape(wcses, 1, []),'filled')
plot(err_av, rT_av, 'k--', 'LineWidth', 1)
xlim([0, 0.5])
ylim([0, 18])
xlabel('mean error', 'fontsize', 14)
ylabel('mean decision time', 'fontsize', 14)
title('       curve(\theta, \alpha)      ', 'fontsize', 14)
set(gcf, 'color', 'w', 'units', 'inches', 'position', [0 0 2.75 2.5])
hold off

if ~exist('figs', 'dir')
       mkdir('figs')
    end

export_fig(append('./figs/', mfilename), '-pdf', '-eps', '-q101');
savefig([pwd '/figs/' mfilename])