clear all; close all; clc; dbstop if error

% please note: due to the stochasticity of generated decision trajectories, there
% may be some varyation in figures produced 

load('3N_curve_data_dense.mat');
N = 3;

% distribution
dist.generate = @generateLogp;
dist.mu = 1/2; dist.si = 1; dist.N = N;
dist.ntrials = 1000000;
dist.maxt = 30;

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

[x1, x2] = meshgrid( linspace(0, 1, opt.nt) );
x1 = opt.rx(1,1) + x1*diff(opt.rx(1,:));
x2 = opt.rx(2,1) + x2*diff(opt.rx(2,:));
thetas = exp(repmat(x1, [1,1,length(reward.wcs)]));
thetas(isnan(r_opt))=nan;
alphas = (repmat(x2, [1,1,length(reward.wcs)]));
alphas(isnan(r_opt))=nan;
 
wcss = [50, 200, 480];
i = 1;
for wcsss = wcss
    thetas_temp = reshape(thetas(:,:,wcsss), 1, []);
    alphas_temp = reshape(alphas(:,:,wcsss), 1, []);
    
    thetas_temp = thetas_temp(~isnan(alphas_temp));
    alphas_temp = alphas_temp(~isnan(alphas_temp));
    
    alphas_temp = alphas_temp(~isnan(thetas_temp));
    thetas_temp = thetas_temp(~isnan(thetas_temp));
    
    %flat, alpha == 0; max theta, or near as poss...
    if (alphas_temp==0); x(2,13,i) = 0; else x(2,[1, 11],i) = [min(abs(alphas_temp), [], 'omitnan'), -min(abs(alphas_temp), [], 'omitnan')]; end
    if (alphas_temp==0); x(1,13,i) = max(thetas_temp(alphas_temp==0), [], 'omitnan'); else; x(1,[13 23],i) = [max(thetas_temp(alphas_temp == x(2,1,i)), [], 'omitnan') max(thetas_temp(alphas_temp == x(2,11,i)), [], 'omitnan')]; end
    
    %flat, alpha == 0; min theta, or near as poss...
    if (alphas_temp==0); x(2,14,i) = 0; else x(2,[2, 12],i) = [min(abs(alphas_temp), [], 'omitnan'), -min(abs(alphas_temp), [], 'omitnan')]; end
    if (alphas_temp==0); x(1,14,i) = min(thetas_temp(alphas_temp==0), [], 'omitnan'); else; x(1,[14 24],i) = [min(thetas_temp(alphas_temp == x(2,2,i)), [], 'omitnan') min(thetas_temp(alphas_temp == x(2,12,i)), [], 'omitnan')]; end
    
    %max theta, max alpha
    x(1,15,i) = max(thetas_temp, [], 'omitnan');
    x(2,15,i) = max(alphas_temp(thetas_temp == x(1,15,i)), [], 'omitnan');
    
    %max theta, min alpha
    x(1,16,i) = max(thetas_temp, [], 'omitnan');
    x(2,16,i) = min(alphas_temp(thetas_temp == x(1,16,i)), [], 'omitnan');
    
    %min theta, max alpha
    x(1,17,i) = min(thetas_temp, [], 'omitnan');
    x(2,17,i) = max(alphas_temp(thetas_temp == x(1,17,i)), [], 'omitnan');
    
    %min theta, min alpha
    x(1,18,i) = min(thetas_temp, [], 'omitnan');
    x(2,18,i) = min(alphas_temp(thetas_temp == x(1,18,i)), [], 'omitnan');
    
    %max alpha, max theta
    x(2,19,i) = max(alphas_temp, [], 'omitnan');
    x(1,19,i) = max(thetas_temp(alphas_temp == x(2,19,i)), [], 'omitnan');
    
    %max alpha, min theta
    x(2,20,i) = max(alphas_temp, [], 'omitnan');
    x(1,20,i) = min(thetas_temp(alphas_temp == x(2,20,i)), [], 'omitnan');
    
    %min alpha, max theta
    x(2,21,i) = min(alphas_temp, [], 'omitnan');
    x(1,21,i) = max(thetas_temp(alphas_temp == x(2,21,i)), [], 'omitnan');
    
    %min alpha, min theta
    x(2,22,i) = min(alphas_temp, [], 'omitnan');
    x(1,22,i) = min(thetas_temp(alphas_temp == x(2,22,i)), [], 'omitnan');
    
    i = i + 1;
end

wcss = [50, 200, 480];
i = 1;
for wcsss = wcss
    thetas_temp = reshape(thetas(:,:,wcsss), 1, []);
    alphas_temp = reshape(alphas(:,:,wcsss), 1, []);
    
    thetas_temp = thetas_temp(~isnan(alphas_temp));
    alphas_temp = alphas_temp(~isnan(alphas_temp));
    
    alphas_temp = alphas_temp(~isnan(thetas_temp));
    thetas_temp = thetas_temp(~isnan(thetas_temp));
    
    x(1,1:12,i) = datasample(thetas_temp, 12);
    
    for j = 1:1:12
        x(2,j,i) = datasample(alphas_temp(x(1,j,i)==(thetas_temp)), 1);
    end
    
    i = i + 1;
end

nsamples = 10000;
 
%% Analyis
 
% generate distribution data
[dist.z,dist.target] = dist.generate(dist);

itrial = randi(dist.ntrials,[nsamples 1]);
 
% curve
scale = @(p, alpha) 1 + (1/nchoosek(size(p, 2), 3)) .* alpha .* sum( prod( reshape(p(:, nchoosek(1:size(p, 2), 3)), [size(p, 1), nchoosek(size(p, 2), 3), 3]), 3), 2);
 
% normalize weights
c = max(wcs, [], 2); ws = wcs./c;
 
% threshold crossing
[e, rt] = deal(nan(nsamples,1));
for isample = 1:nsamples
    for combos = 1:size(x, 2)
        for wces = 1:size(x,3)
            sc = scale( exp( dist.z( 2:end, :, itrial(isample) ) ), x(2,combos,wces) );
            sc = [sc(1:end-1); 0]; % solve crossing problem
            isdec = exp(dist.z(2:end, :, itrial(isample))) > x(1,combos,wces)*sc.*ones(1,N);
            rt(isample,combos,wces) = find( sum( isdec ,2), 1);
            [~, n] = max( dist.z(rt(isample,combos,wces)+1, :, itrial(isample)) );
            sc = scale( exp( dist.z( find( sum( isdec ,2), 1), :, itrial(isample) ) ), x(2,combos,wces));
            bound(isample,combos,wces) = x(1,combos,wces)*sc;
            e(isample,combos,wces) = n~=dist.target(itrial(isample));
        end
    end
end
for combos = 1:size(x, 2)
    for wces = 1:size(x,3)
        rnan = rt(:,combos,wces)>dist.maxt; e(rnan,combos,wces) = nan; rt(rnan,combos,wces) = nan; bound(rnan,combos,wces) = nan;
    end
end

for wces = 1:size(x,3)
   for combos = 1:size(x,2)

       [rt_b_avg{combos,wces},ic1,rnk]=unique(rt(:,combos,wces));
       
       for r_t = 1:max(rnk)
           a = [mean(bound(rnk==r_t,combos,wces)); max(bound(rnk==r_t,combos,wces)); min(bound(rnk==r_t,combos,wces))];
           b = log(a./(1-a));
           [a; b];
           
            bound_avg{combos,wces}(r_t)=mean(bound(rnk==r_t,combos,wces));
            bound_var{combos,wces}(r_t)=var(bound(rnk==r_t,combos,wces));
            bound_max{combos,wces}(r_t)=max(bound(rnk==r_t,combos,wces));
            bound_min{combos,wces}(r_t)=min(bound(rnk==r_t,combos,wces));
            
            log_bound = log(bound(rnk==r_t,combos,wces)./(1-bound(rnk==r_t,combos,wces)));
            
            %average of transformed
            bound_avg_logd{combos,wces}(r_t)=mean(log_bound);
            bound_var_logd{combos,wces}(r_t)=var(log_bound);
            bound_max_logd{combos,wces}(r_t)=max(log_bound);
            bound_min_logd{combos,wces}(r_t)=min(log_bound);
            
            %average, THEN transformed
            bound_avg_logd{combos,wces}(r_t) = log((bound_avg{combos,wces}(r_t))./(1-(bound_avg{combos,wces}(r_t))));
            
       end
   end
end


for wces = 1:size(x,3)
   for combos = 1:size(x,2)
       scatter(rt_b_avg{combos,wces}, bound_avg_logd{combos,wces}, 20, 'filled')
       title("\theta = " + x(1,combos,wces) + ", \alpha = " + x(2, combos, wces) + "\n Combos = " + combos + ", wces = " + wces)
       xlim([0 30])
       ylim([-0.5 6.5])
   end
end

%%

N = 3;
C = linspecer(N);
C = [30, 136, 229; 216, 27, 96; 255, 193, 7]./255;

i=1;
%flat
combos_i = [13, 14, 13, 14, 12, 13];
wces = [1, 1, 2, 2, 3, 3];
figure
for combos = combos_i
    scatter(rt_b_avg{combos,wces(i)}(2:end), bound_avg_logd{combos,wces(i)}(2:end), 10, C(wces(i),:), 'filled')
    hold on
    
    gprMdl_mn = fitrgp(rt_b_avg{combos,wces(i)}(2:end), real(bound_avg_logd{combos,wces(i)}(2:end)),'Basis','pureQuadratic','FitMethod','exact',...
        'PredictMethod','exact','SigmaLowerBound',0.251);
    ypred = resubPredict(gprMdl_mn);
    X = rt_b_avg{combos,wces(i)}(~isnan(real(bound_avg_logd{combos,wces(i)})));
    plot(X(2:end), ypred, 'Color', C(wces(i),:), 'LineWidth', 2)
    xlim([0 30])
    ylim([-0.5 6.5])
    hold on
    rt_f=[rt_f; rt_f];
    bnds_f=[bnds_f, bnds_f];
i=i+1;
end

xlabel('reaction time ({\it rT})')
hf = gcf;
ha = gca;
set(hf, 'color', 'w', 'units', 'inches', 'position', [0 0 4 4])
set(ha, 'color', 'w', 'units', 'inches', 'position', [1 1 2.75 2.5])
ha.XAxis.TickValues = [0 10 20 30];
ha.YAxis.TickValues = [0 2 4 6];

hyl = ylabel({'{\bf \fontsize{20}Static}','average transformed threshold'}); hyl.Position = hyl.Position+[0.3 0 0];

if ~exist('figs', 'dir')
       mkdir('figs')
    end

export_fig(append('./figs/fig8g'), '-pdf', '-eps', '-q101');
savefig([pwd '/figs/fig8g'])

%%
i=1;
%decreasing
combos_i = [7, 10, 3, 15, 15];
wces = [1, 1, 2, 2, 3];
figure
for combos = combos_i
    scatter(rt_b_avg{combos,wces(i)}(2:end), bound_avg_logd{combos,wces(i)}(2:end), 10, C(wces(i),:), 'filled')
    hold on
    
    gprMdl_mn = fitrgp(rt_b_avg{combos,wces(i)}(2:end), real(bound_avg_logd{combos,wces(i)}(2:end)),'Basis','pureQuadratic','FitMethod','exact',...
        'PredictMethod','exact','SigmaLowerBound',0.251);
    ypred = resubPredict(gprMdl_mn);
    X = rt_b_avg{combos,wces(i)}(~isnan(real(bound_avg_logd{combos,wces(i)})));
    plot(X(2:end), ypred, 'Color', C(wces(i),:), 'LineWidth', 2)
    xlim([0 30])
    ylim([-0.5 6.5])
    hold on
    rt_d=[rt_d; rt_d];
    bnds_d=[bnds_d, bnds_d];
i=i+1;
end
hf = gcf;
ha = gca;
set(hf, 'color', 'w', 'units', 'inches', 'position', [0 0 4 4])
set(ha, 'color', 'w', 'units', 'inches', 'position', [1 1 2.75 2.5])
ha.XAxis.TickValues = [0 10 20 30];
ha.YAxis.TickValues = [0 2 4 6];

hyl = ylabel({'{\bf \fontsize{20}Decreasing}','average transformed threshold'}); hyl.Position = hyl.Position+[0.3 0 0];

if ~exist('figs', 'dir')
       mkdir('figs')
    end

export_fig(append('./figs/fig8d'), '-pdf', '-eps', '-q101');
savefig([pwd '/figs/fig8d'])

%%
i=1;
j=1;
%increasing
C = [30, 136, 229; 216, 27, 96; 255, 193, 7]./255;
combos_i = [22, 2, 3, 1, 9, 10, 5];
wces = [1, 1, 1, 2, 2, 2, 2];
figure
for combos = combos_i
    scatter(rt_b_avg{combos,wces(i)}(2:end), bound_avg_logd{combos,wces(i)}(2:end), 10, C(wces(i),:), 'filled')
    hold on
    
    gprMdl_mn = fitrgp(rt_b_avg{combos,wces(i)}(2:end), real(bound_avg_logd{combos,wces(i)}(2:end)),'Basis','pureQuadratic','FitMethod','exact',...
        'PredictMethod','exact','SigmaLowerBound',0.251);
    ypred = resubPredict(gprMdl_mn);
    X = rt_b_avg{combos,wces(i)}(~isnan(real(bound_avg_logd{combos,wces(i)})));
    plot(X(2:end), ypred, 'Color', C(wces(i),:), 'LineWidth', 2)
    xlim([0 30])
    ylim([-0.5 6.5])
    hold on
    rt_i=[rt_i; rt_i];
    bnds_i=[bnds_i, bnds_i];
i=i+1;
end
hf = gcf;
ha = gca;
set(hf, 'color', 'w', 'units', 'inches', 'position', [0 0 4 4])
set(ha, 'color', 'w', 'units', 'inches', 'position', [1 1 2.75 2.5])
title(ha, 'curve(\theta, \alpha)', 'FontSize', 20);
ha.XAxis.TickValues = [0 10 20 30];
ha.YAxis.TickValues = [0 2 4 6];

hyl = ylabel({'{\bf \fontsize{20}Increasing}','average transformed threshold'}); hyl.Position = hyl.Position+[0.3 0 0];

if ~exist('figs', 'dir')
       mkdir('figs')
    end

export_fig(append('./figs/fig8a'), '-pdf', '-eps', '-q101');
savefig([pwd '/figs/fig8a'])