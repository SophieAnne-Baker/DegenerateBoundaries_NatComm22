%% combine createFitsWcsn_dense

close all; clear all; clf

%% load data/fits
[datx1, daty1, fitresult1, fts1, optss1, thetas1, alphas1, r_opt1] = createFitsWcs1_dense(0,0);

[datx2, daty2, fitresult2, fts2, optss2, thetas2, alphas2, r_opt2] = createFitsWcs2_dense(0,0);

[datx3, daty3, fitresult3, fts3, optss3, thetas3, alphas3, r_opt3] = createFitsWcs3_dense(0,0);

%% Plot fitresult(s) data
cmap = parula(5);
sz = 10;
h_adj =[-0.07, 0.08];

figure;
set(gcf, 'color', 'w', 'units', 'inches', 'position', [0 0 12 11])
subplot(4,2,2); hold on
for i=1:5
    h1{i} = plot( fitresult1{i}, datx1{i}, daty1{i}, 'predobs' );
    h1{i}(1).Color = cmap(i,:); %data
    h1{i}(2).Color = cmap(i,:); %linefit
    h1{i}(2).LineWidth = 2.5;
    h1{i}(3).Color = cmap(i,:); %lower conf bound
    h1{i}(4).Color = cmap(i,:); %upper conf bound
    p1 = patch([h1{i}(3).XData fliplr(h1{i}(4).XData)], [h1{i}(3).YData fliplr(h1{i}(4).YData)], cmap(i,:));
    p1.FaceAlpha = 0.1;
    peak_lower(i) = nanmax(h1{i}(3).YData);
    peak_upper(i) = nanmax(h1{i}(4).YData);
end
legend('off')
peak_X_min = nanmin([h1{1}(3).XData, h1{2}(3).XData, h1{3}(3).XData, h1{4}(3).XData, h1{5}(3).XData])-0.01;
peak_X_max = nanmax([h1{1}(4).XData, h1{2}(4).XData, h1{3}(4).XData, h1{4}(4).XData, h1{5}(4).XData])+0.009;
plot(linspace(peak_X_min, peak_X_max, 100), max(peak_lower)*ones(1,100), '--r', 'LineWidth', 2);
plot(linspace(peak_X_min, peak_X_max, 100), min(peak_upper)*ones(1,100), '--r', 'LineWidth', 2);
xlim([0.9, 0.99])
title('')
xlabel('\theta','fontsize',15)
ylabel('average reward,{\it r}','fontsize',15)
title({['Sections'],[' ']},'fontsize',20);

subplot(4,2,4); hold on
for i=1:5
    h2{i} = plot( fitresult2{i}, datx2{i}, daty2{i}, 'predobs' );
    h2{i}(1).Color = cmap(i,:);
    h2{i}(2).Color = cmap(i,:);
    h2{i}(2).LineWidth = 2.5;
    h2{i}(3).Color = cmap(i,:);
    h2{i}(4).Color = cmap(i,:);
    p2 = patch([h2{i}(3).XData fliplr(h2{i}(4).XData)], [h2{i}(3).YData fliplr(h2{i}(4).YData)], cmap(i,:));
    p2.FaceAlpha = 0.1;
    peak_lower(i) = nanmax(h2{i}(3).YData);
    peak_upper(i) = nanmax(h2{i}(4).YData);
end
legend('off')
peak_X_min = nanmin([h2{1}(3).XData, h2{2}(3).XData, h2{3}(3).XData, h2{4}(3).XData, h2{5}(3).XData])-0.01;
peak_X_max = nanmax([h2{1}(4).XData, h2{2}(4).XData, h2{3}(4).XData, h2{4}(4).XData, h2{5}(4).XData])+0.009;
plot(linspace(peak_X_min, peak_X_max, 100), max(peak_lower)*ones(1,100), '--r', 'LineWidth', 2);
plot(linspace(peak_X_min, peak_X_max, 100), min(peak_upper)*ones(1,100), '--r', 'LineWidth', 2);
xlim([0.57, 0.85])
xlabel('\theta','fontsize',15)
ylabel('average reward,{\it r}','fontsize',15)

clear peak_lower peak_upper peak_X_min peak_X_max
subplot(4,2,6); hold on
for i=2:5
    h3{i} = plot( fitresult3{i}, datx3{i}, daty3{i}, 'predobs' );
    h3{i}(1).Color = cmap(i,:);
    h3{i}(2).Color = cmap(i,:);
    h3{i}(2).LineWidth = 2.5;
    h3{i}(3).Color = cmap(i,:);
    h3{i}(4).Color = cmap(i,:);
    p3 = patch([h3{i}(3).XData fliplr(h3{i}(4).XData)], [h3{i}(3).YData fliplr(h3{i}(4).YData)], cmap(i,:));
    p3.FaceAlpha = 0.1;
    peak_lower(i-1) = nanmax(h3{i}(3).YData);
    peak_upper(i-1) = nanmax(h3{i}(4).YData);
end
legend('off')
peak_X_min = nanmin([h3{2}(3).XData, h3{3}(3).XData, h3{4}(3).XData, h3{5}(3).XData])-0.01;
peak_X_max = nanmax([h3{2}(4).XData, h3{3}(4).XData, h3{4}(4).XData, h3{5}(4).XData])+0.009;
plot(linspace(peak_X_min, peak_X_max, 100), max(peak_lower)*ones(1,100), '--r', 'LineWidth', 2);
plot(linspace(peak_X_min, peak_X_max, 100), min(peak_upper)*ones(1,100), '--r', 'LineWidth', 2);
xlim([0.3, 0.7])
xlabel('\theta','fontsize',15)
ylabel('average reward,{\it r}','fontsize',15)

%% plot landscapes

clear all

sz = 10;
h_adj =[-0.07, 0.08];

load('3N_curve_data_banded.mat');

%extract optimal params
for wcs = 1:1:length(reward.wcs)
    r_temp = r(:,:,wcs);
    r_temp(r_temp<-1) = nan;
    r_std_temp = r_std(:,:,wcs);
    r_std_temp(r_std_temp<-1) = nan;
    del = 0.03.*r_std_temp;
    r_temp(r_temp < max(max(r_temp)) - del) = nan;
    r_opt(:,:,wcs) = r_temp;
end
% get error and RTs for each optimal value
err = repmat(squeeze(v(:,:,1)), [1,1,length(reward.wcs)]);
rT = repmat(squeeze(v(:,:,2)), [1,1,length(reward.wcs)]);
err(isnan(r_opt))=nan;
rT(isnan(r_opt))=nan;

err_std = repmat(squeeze(v_std(:,:,1)), [1,1,length(reward.wcs)]);
rT_std = repmat(squeeze(v_std(:,:,2)), [1,1,length(reward.wcs)]);
err_std(isnan(r_opt))=nan;
rT_std(isnan(r_opt))=nan;
clear r_temp

thetas = exp(repmat(x1, [1,1,length(reward.wcs)]));
thetas(isnan(r_opt))=nan;
alphas = (repmat(x2, [1,1,length(reward.wcs)]));
alphas(isnan(r_opt))=nan;

cw_vals{1} = ['\bf', '0.001'];
cw_vals{2} = ['\bf', '0.04'];
cw_vals{3} = ['\bf', '0.1'];
cw_vals = [0.001, 0.04, 0.1];
for j = 1:3
    sp{j} = subplot(4,2,2*j-1); hold on
    r(r<-1) = nan;
    mesh(exp(x1), x2, r(:,:,j));
    axis([0.4 1 -20 20 -1 -0]);
    set(gca,'xtick',0.4:0.1:1,'ytick',-20:10:20,'ztick',-1:0.1:0)
    set(gca,'fontsize',8)
    color_map=(bone(512));
    color_map = color_map(1:470,:);
    colormap(color_map)
    hc = colorbar; 
    hyl = ylabel(hc,'reward,{\it r}','fontsize',15); 
    hx = xlabel('\Theta','fontsize',15,'rotation',20);
    hy = ylabel('\alpha','fontsize',15,'rotation',-30);
    zlabel({append('{\fontsize{20} \bf {\it c/W} = ',num2str(cw_vals(j)),'}');'';'';'{\fontsize{15} average reward,{\it r}}'});
    if(j==1); title({['Reward Landscapes'],[' ']},'fontsize',20); end

    view([-35 20])
    cmap = parula(5);
    for i=1:5
        thscat = thetas(i,:,j); thscat = thscat(:);
        alphscat = ones(size(thscat));
        z = ones(size(thscat));
        scatter3(thscat, x2(i,1).*alphscat, -z, sz, cmap(i,:), 'filled')
        scatter3(thscat, x2(i,1).*alphscat+1, -z, sz, cmap(i,:), 'filled');
        scatter3(thscat, x2(i,1).*alphscat+0.5, -z, sz, cmap(i,:), 'filled');
        scatter3(thscat, x2(i,1).*alphscat, r_opt(i,:,j), sz, cmap(i,:), 'filled');
        scatter3(thscat, x2(i,1).*alphscat+1, r_opt(i,:,j), sz, cmap(i,:), 'filled');
        scatter3(thscat, x2(i,1).*alphscat+0.5, r_opt(i,:,j), sz, cmap(i,:), 'filled');
        
        
    end

    scatter3(thscat, x2(5,1).*alphscat-0.5, r_opt(5,:,j), sz, cmap(5,:), 'filled');
    scatter3(thscat, x2(5,1).*alphscat-1, r_opt(5,:,j), sz, cmap(5,:), 'filled');
    scatter3(thscat, x2(5,1).*alphscat+1.5, r_opt(5,:,j), sz, cmap(5,:), 'filled');
    scatter3(thscat, x2(5,1).*alphscat+2, r_opt(5,:,j), sz, cmap(5,:), 'filled');
    scatter3(thscat, x2(5,1).*alphscat-0.5, -z, sz, cmap(5,:), 'filled');
    scatter3(thscat, x2(5,1).*alphscat-1, -z, sz, cmap(5,:), 'filled');
    scatter3(thscat, x2(5,1).*alphscat+1.5, -z, sz, cmap(5,:), 'filled');
    scatter3(thscat, x2(5,1).*alphscat+2, -z, sz, cmap(5,:), 'filled');
    
    %smooth projection
    z_mtrx = ones(size(x2));
    mesh(exp(x1), x2, -z_mtrx, r(:,:,j));

end



%% extract optimal params, colour optima in bounds
for wcs = 1:1:length(reward.wcs)
    r_temp = r(:,:,wcs);
    r_temp(r_temp<-1) = nan;
    r_std_temp = r_std(:,:,wcs);
    r_std_temp(r_std_temp<-1) = nan;
    del = 0.01.*r_std_temp;
    r_temp(r_temp < max(max(r_temp)) - del) = nan;
    r_opt(:,:,wcs) = r_temp;
end
% get error and RTs for each optimal value
err = repmat(squeeze(v(:,:,1)), [1,1,length(reward.wcs)]);
rT = repmat(squeeze(v(:,:,2)), [1,1,length(reward.wcs)]);
err(isnan(r_opt))=nan;
rT(isnan(r_opt))=nan;

err_std = repmat(squeeze(v_std(:,:,1)), [1,1,length(reward.wcs)]);
rT_std = repmat(squeeze(v_std(:,:,2)), [1,1,length(reward.wcs)]);
err_std(isnan(r_opt))=nan;
rT_std(isnan(r_opt))=nan;
clear r_temp

thetas = exp(repmat(x1, [1,1,length(reward.wcs)]));
thetas(isnan(r_opt))=nan;
alphas = (repmat(x2, [1,1,length(reward.wcs)]));
alphas(isnan(r_opt))=nan;

a = subplot(4,2,2); hold on
j=1; for i=1:5
thscat = thetas(i,:,j); thscat = thscat(:);
alphscat = ones(size(thscat));
z = ones(size(thscat));
scatter(thscat, r(i,:,j), sz,'k', 'filled')
end
b = subplot(4,2,4); hold on
j=2; for i=1:5
thscat = thetas(i,:,j); thscat = thscat(:);
alphscat = ones(size(thscat));
z = ones(size(thscat));
scatter(thscat, r(i,:,j), sz,'k', 'filled')
end
c = subplot(4,2,6); hold on
j=3; for i=1:5
thscat = thetas(i,:,j); thscat = thscat(:);
alphscat = ones(size(thscat));
z = ones(size(thscat));
scatter(thscat, r(i,:,j), sz,'k', 'filled')
end


d = subplot(4,2,1); hold on
j=1; for i=1:5
thscat = thetas(i,:,j); thscat = thscat(:);
alphscat = ones(size(thscat));
z = ones(size(thscat))-0.01;
scatter3(thscat, x2(i,1).*alphscat, -z, sz,'k', 'filled')
scatter3(thscat, x2(i,1).*alphscat, -z, sz, 'k', 'filled')
scatter3(thscat, x2(i,1).*alphscat+1, -z, sz, 'k', 'filled');
scatter3(thscat, x2(i,1).*alphscat+0.5, -z, sz, 'k', 'filled');

scatter3(thscat, x2(i,1).*alphscat, r(i,:,j), sz,'k', 'filled')
scatter3(thscat, x2(i,1).*alphscat, r(i,:,j), sz, 'k', 'filled')
scatter3(thscat, x2(i,1).*alphscat+1, r(i,:,j), sz, 'k', 'filled');
scatter3(thscat, x2(i,1).*alphscat+0.5, r(i,:,j), sz, 'k', 'filled');
end
scatter3(thscat, x2(5,1).*alphscat-0.5, -z, sz, 'k', 'filled');
scatter3(thscat, x2(5,1).*alphscat-1, -z, sz, 'k', 'filled');
scatter3(thscat, x2(5,1).*alphscat+1.5, -z, sz, 'k', 'filled');
scatter3(thscat, x2(5,1).*alphscat+2, -z, sz, 'k', 'filled');

scatter3(thscat, x2(5,1).*alphscat-0.5, r(5,:,j), sz, 'k', 'filled');
scatter3(thscat, x2(5,1).*alphscat-1, r(5,:,j), sz, 'k', 'filled');
scatter3(thscat, x2(5,1).*alphscat+1.5, r(5,:,j), sz, 'k', 'filled');
scatter3(thscat, x2(5,1).*alphscat+2, r(5,:,j), sz, 'k', 'filled');

e = subplot(4,2,3); hold on
j=2; for i=1:5
thscat = thetas(i,:,j); thscat = thscat(:);
alphscat = ones(size(thscat));
z = ones(size(thscat))-0.01;
scatter3(thscat, x2(i,1).*alphscat, -z, sz,'k', 'filled')
scatter3(thscat, x2(i,1).*alphscat, -z, sz, 'k', 'filled')
scatter3(thscat, x2(i,1).*alphscat+1, -z, sz, 'k', 'filled');
scatter3(thscat, x2(i,1).*alphscat+0.5, -z, sz, 'k', 'filled');

scatter3(thscat, x2(i,1).*alphscat, r(i,:,j), sz,'k', 'filled')
scatter3(thscat, x2(i,1).*alphscat, r(i,:,j), sz, 'k', 'filled')
scatter3(thscat, x2(i,1).*alphscat+1, r(i,:,j), sz, 'k', 'filled');
scatter3(thscat, x2(i,1).*alphscat+0.5, r(i,:,j), sz, 'k', 'filled');
end
scatter3(thscat, x2(5,1).*alphscat-0.5, -z, sz, 'k', 'filled');
scatter3(thscat, x2(5,1).*alphscat-1, -z, sz, 'k', 'filled');
scatter3(thscat, x2(5,1).*alphscat+1.5, -z, sz, 'k', 'filled');
scatter3(thscat, x2(5,1).*alphscat+2, -z, sz, 'k', 'filled');

scatter3(thscat, x2(5,1).*alphscat-0.5, r(5,:,j), sz, 'k', 'filled');
scatter3(thscat, x2(5,1).*alphscat-1, r(5,:,j), sz, 'k', 'filled');
scatter3(thscat, x2(5,1).*alphscat+1.5, r(5,:,j), sz, 'k', 'filled');
scatter3(thscat, x2(5,1).*alphscat+2, r(5,:,j), sz, 'k', 'filled');

f = subplot(4,2,5); hold on
j=3; for i=1:5
thscat = thetas(i,:,j); thscat = thscat(:);
alphscat = ones(size(thscat));
z = ones(size(thscat))-0.01;
scatter3(thscat, x2(i,1).*alphscat, -z, sz,'k', 'filled')
scatter3(thscat, x2(i,1).*alphscat, -z, sz, 'k', 'filled')
scatter3(thscat, x2(i,1).*alphscat+1, -z, sz, 'k', 'filled');
scatter3(thscat, x2(i,1).*alphscat+0.5, -z, sz, 'k', 'filled');

scatter3(thscat, x2(i,1).*alphscat, r(i,:,j), sz,'k', 'filled')
scatter3(thscat, x2(i,1).*alphscat, r(i,:,j), sz, 'k', 'filled')
scatter3(thscat, x2(i,1).*alphscat+1, r(i,:,j), sz, 'k', 'filled');
scatter3(thscat, x2(i,1).*alphscat+0.5, r(i,:,j), sz, 'k', 'filled');
end
scatter3(thscat, x2(5,1).*alphscat-0.5, -z, sz, 'k', 'filled');
scatter3(thscat, x2(5,1).*alphscat-1, -z, sz, 'k', 'filled');
scatter3(thscat, x2(5,1).*alphscat+1.5, -z, sz, 'k', 'filled');
scatter3(thscat, x2(5,1).*alphscat+2, -z, sz, 'k', 'filled');

scatter3(thscat, x2(5,1).*alphscat-0.5, r(5,:,j), sz, 'k', 'filled');
scatter3(thscat, x2(5,1).*alphscat-1, r(5,:,j), sz, 'k', 'filled');
scatter3(thscat, x2(5,1).*alphscat+1.5, r(5,:,j), sz, 'k', 'filled');
scatter3(thscat, x2(5,1).*alphscat+2, r(5,:,j), sz, 'k', 'filled');

a.Position([2,4]) = a.Position([2,4])+h_adj;
b.Position([2,4]) = b.Position([2,4])+[2*h_adj(1), h_adj(2)];
c.Position([2,4]) = c.Position([2,4])+[3*h_adj(1), h_adj(2)];

d.Position([2,4]) = d.Position([2,4])+h_adj;
e.Position([2,4]) = e.Position([2,4])+[2*h_adj(1), h_adj(2)];
f.Position([2,4]) = f.Position([2,4])+[3*h_adj(1), h_adj(2)];

subplot(4,1,4).Position(2) = 4*h_adj(1);
C = ([parula(5); 0 0 0; 1 0 0]);
lines = {'-','-','-','-','-','.','--'};
l = legend;
for i = 1:5
p{i} = plot(nan(2,1),nan(2,1), 'LineStyle', lines{i}, 'color', C(i,:), 'LineWidth', 2); hold on;
end
p{6} = scatter(nan, nan, 10,'k', 'filled'); hold on;
i=7;
p{7} = plot(nan(2,1),nan(2,1), '--r', 'LineWidth', 2);
l=legend([p{1}, p{2}, p{3}, p{4}, p{5}, p{6}, p{7}], {'\alpha = -20', '\alpha = -10', '\alpha = 0', '\alpha = 10', '\alpha = 20', 'extracted maxima', 'shared confidence range'}, 'Orientation', 'horizontal', 'FontSize', 12, 'Position', [0.0891 0.045 0.8067 0.0257]);

%%

if ~exist('figs', 'dir')
       mkdir('figs')
    end

export_fig(append('./figs/', mfilename), '-pdf', '-eps', '-q101');
savefig([pwd '/figs/' mfilename])