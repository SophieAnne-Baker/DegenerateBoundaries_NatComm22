clear all; 
close all; clc; 
dbstop if error

% please note: due to the stochasticity of generated decision trajectories, there
% may be some varyation in figures produced 

load('3N_power_data.mat');
N = 3;

% distribution
dist.generate = @generateLogp;
dist.mu = 1/2; dist.si = 1; dist.N = N;
dist.ntrials = 100000;
dist.maxt = 30;

%extract optimal params
for wcs = 1:1:length(reward.wcs)
    r_temp = r(:,:,:,wcs);
    r_temp(r_temp<-1) = nan;
    r_std_temp = r_std(:,:,:,wcs);
    r_std_temp(r_std_temp<-1) = nan;
    del = 0.03.*r_std_temp;
    r_temp(r_temp < max(max(r_temp)) - del) = nan;
    r_opt(:,:,:,wcs) = r_temp;
end
% get error and RTs for each optimal value
err = repmat(squeeze(v(:,:,:,1)), [1,1,length(reward.wcs)]);
rT = repmat(squeeze(v(:,:,:,2)), [1,1,length(reward.wcs)]);
err(isnan(r_opt))=nan;
rT(isnan(r_opt))=nan;
clear r_temp

[x1, x2, x3] = meshgrid( linspace(0, 1, opt.nt) );
x1 = opt.rx(1,1) + x1*diff(opt.rx(1,:));
x2 = opt.rx(2,1) + x2*diff(opt.rx(2,:));
x3 = opt.rx(3,1) + x3*diff(opt.rx(3,:));
thetas = exp(repmat(x1, [1,1,1,length(reward.wcs)]));
thetas(isnan(r_opt))=nan;
alphas = (repmat(x2, [1,1,1,length(reward.wcs)]));
alphas(isnan(r_opt))=nan;
betas = (repmat(x3, [1,1,1,length(reward.wcs)]));
betas(isnan(r_opt))=nan;

wcss = [50, 200, 480];
i = 1;
for wcsss = wcss
    thetas_temp = reshape(thetas(:,:,:,wcsss), 1, []);
    alphas_temp = reshape(alphas(:,:,:,wcsss), 1, []);
    betas_temp = reshape(betas(:,:,:,wcsss), 1, []);
    
    thetas_temp = thetas_temp(~isnan(betas_temp));
    alphas_temp = alphas_temp(~isnan(betas_temp));
    betas_temp = betas_temp(~isnan(betas_temp));
    
    betas_temp = betas_temp(~isnan(alphas_temp));
    thetas_temp = thetas_temp(~isnan(alphas_temp));
    alphas_temp = alphas_temp(~isnan(alphas_temp));
    
    alphas_temp = alphas_temp(~isnan(thetas_temp));
    betas_temp = betas_temp(~isnan(thetas_temp));
    thetas_temp = thetas_temp(~isnan(thetas_temp));
    
    %flat, alpha == 0; max theta, or near as poss... 
    if (alphas_temp==0)
        x(2,1,i) = 0;
        x(1,1,i) = thetas_temp(alphas_temp==0);
        x(3,1,i) = 0;
    else
        % min alpha 0+, min theta, min beta
        x(2,1,i) = min(abs(alphas_temp), [], 'omitnan');
        x(1,1,i) = min(thetas_temp(alphas_temp==x(2,1,i)), [], 'omitnan');
        x(3,1,i) = min(betas_temp(alphas_temp==x(2,1,i) & thetas_temp == x(1,1,i)), [], 'omitnan' );
        
        % min alpha 0+, min theta, max beta
        x(2,50,i) = min(abs(alphas_temp), [], 'omitnan');
        x(1,50,i) = min(thetas_temp(alphas_temp==x(2,50,i)), [], 'omitnan');
        x(3,50,i) = max(betas_temp(alphas_temp==x(2,50,i) & thetas_temp == x(1,50,i)), [], 'omitnan' );
        
        % min alpha 0+, max theta, min beta
        x(2,51,i) = min(abs(alphas_temp), [], 'omitnan');
        x(1,51,i) = max(thetas_temp(alphas_temp==x(2,51,i)), [], 'omitnan');
        x(3,51,i) = min(betas_temp(alphas_temp==x(2,51,i) & thetas_temp == x(1,51,i)), [], 'omitnan' );
        
        % min alpha 0+, max theta, max beta
        x(2,52,i) = min(abs(alphas_temp), [], 'omitnan');
        x(1,52,i) = max(thetas_temp(alphas_temp==x(2,52,i)), [], 'omitnan');
        x(3,52,i) = max(betas_temp(alphas_temp==x(2,52,i) & thetas_temp == x(1,52,i)), [], 'omitnan' );
        
        % min alpha 0+, min beta, min theta
        x(2,53,i) = min(abs(alphas_temp), [], 'omitnan');
        x(3,53,i) = min(betas_temp(alphas_temp==x(2,53,i)), [], 'omitnan');
        x(1,53,i) = min(thetas_temp(alphas_temp==x(2,53,i) & betas_temp == x(3,53,i)), [], 'omitnan' );
        
        % min alpha 0+, min beta, max theta
        x(2,54,i) = min(abs(alphas_temp), [], 'omitnan');
        x(3,54,i) = min(betas_temp(alphas_temp==x(2,54,i)), [], 'omitnan');
        x(1,54,i) = max(thetas_temp(alphas_temp==x(2,54,i) & betas_temp == x(3,54,i)), [], 'omitnan' );
        
        % min alpha 0+, max beta, min theta
        x(2,55,i) = min(abs(alphas_temp), [], 'omitnan');
        x(3,55,i) = max(betas_temp(alphas_temp==x(2,55,i)), [], 'omitnan');
        x(1,55,i) = min(thetas_temp(alphas_temp==x(2,55,i) & betas_temp == x(3,55,i)), [], 'omitnan' );
        
        % min alpha 0+, max beta, max theta
        x(2,56,i) = min(abs(alphas_temp), [], 'omitnan');
        x(3,56,i) = max(betas_temp(alphas_temp==x(2,56,i)), [], 'omitnan');
        x(1,56,i) = max(thetas_temp(alphas_temp==x(2,56,i) & betas_temp == x(3,56,i)), [], 'omitnan' );
        
        % min alpha 0-
    
        % min alpha 0-, min theta, min beta
        x(2,57,i) = -min(abs(alphas_temp), [], 'omitnan');
        x(1,57,i) = min(thetas_temp(alphas_temp==x(2,57,i)), [], 'omitnan');
        x(3,57,i) = min(betas_temp(alphas_temp==x(2,57,i) & thetas_temp == x(1,57,i)), [], 'omitnan' );
        
        % min alpha 0-, min theta, max beta
        x(2,58,i) = -min(abs(alphas_temp), [], 'omitnan');
        x(1,58,i) = min(thetas_temp(alphas_temp==x(2,58,i)), [], 'omitnan');
        x(3,58,i) = max(betas_temp(alphas_temp==x(2,58,i) & thetas_temp == x(1,58,i)), [], 'omitnan' );
        
        % min alpha 0-, max theta, min beta
        x(2,59,i) = -min(abs(alphas_temp), [], 'omitnan');
        x(1,59,i) = max(thetas_temp(alphas_temp==x(2,59,i)), [], 'omitnan');
        x(3,59,i) = min(betas_temp(alphas_temp==x(2,59,i) & thetas_temp == x(1,59,i)), [], 'omitnan' );
        
        % min alpha 0-, max theta, max beta
        x(2,60,i) = -min(abs(alphas_temp), [], 'omitnan');
        x(1,60,i) = max(thetas_temp(alphas_temp==x(2,60,i)), [], 'omitnan');
        x(3,60,i) = max(betas_temp(alphas_temp==x(2,60,i) & thetas_temp == x(1,60,i)), [], 'omitnan' );
        
        % min alpha 0-, min beta, min theta
        x(2,61,i) = -min(abs(alphas_temp), [], 'omitnan');
        x(3,61,i) = min(betas_temp(alphas_temp==x(2,61,i)), [], 'omitnan');
        x(1,61,i) = min(thetas_temp(alphas_temp==x(2,61,i) & betas_temp == x(3,61,i)), [], 'omitnan' );
        
        % min alpha 0-, min beta, max theta
        x(2,62,i) = -min(abs(alphas_temp), [], 'omitnan');
        x(3,62,i) = min(betas_temp(alphas_temp==x(2,62,i)), [], 'omitnan');
        x(1,62,i) = max(thetas_temp(alphas_temp==x(2,62,i) & betas_temp == x(3,62,i)), [], 'omitnan' );
        
        % min alpha 0-, max beta, min theta
        x(2,63,i) = -min(abs(alphas_temp), [], 'omitnan');
        x(3,63,i) = max(betas_temp(alphas_temp==x(2,63,i)), [], 'omitnan');
        x(1,63,i) = min(thetas_temp(alphas_temp==x(2,63,i) & betas_temp == x(3,63,i)), [], 'omitnan' );
        
        % min alpha 0-, max beta, max theta
        x(2,64,i) = -min(abs(alphas_temp), [], 'omitnan');
        x(3,64,i) = max(betas_temp(alphas_temp==x(2,64,i)), [], 'omitnan');
        x(1,64,i) = max(thetas_temp(alphas_temp==x(2,64,i) & betas_temp == x(3,64,i)), [], 'omitnan' );
        
    end
    
    % theta lead %
    
    %max theta, max alpha, max beta
    x(1,2,i) = max(thetas_temp, [], 'omitnan');
    x(2,2,i) = max(alphas_temp(thetas_temp == x(1,2,i)), [], 'omitnan');
    x(3,2,i) = max(betas_temp(thetas_temp == x(1,2,i) & alphas_temp == x(2,2,i)), [], 'omitnan');
    
    %max theta, max alpha, min beta
    x(1,3,i) = max(thetas_temp, [], 'omitnan');
    x(2,3,i) = max(alphas_temp(thetas_temp == x(1,3,i)), [], 'omitnan');
    x(3,3,i) = min(betas_temp(thetas_temp == x(1,3,i) & alphas_temp == x(2,3,i)), [], 'omitnan');
    
    %max theta, min alpha, max beta
    x(1,4,i) = max(thetas_temp, [], 'omitnan');
    x(2,4,i) = min(alphas_temp(thetas_temp == x(1,4,i)), [], 'omitnan');
    x(3,4,i) = max(betas_temp(thetas_temp == x(1,4,i) & alphas_temp == x(2,4,i)), [], 'omitnan');
    
    %max theta, min alpha, min beta
    x(1,5,i) = max(thetas_temp, [], 'omitnan');
    x(2,5,i) = min(alphas_temp(thetas_temp == x(1,5,i)), [], 'omitnan');
    x(3,5,i) = min(betas_temp(thetas_temp == x(1,5,i) & alphas_temp == x(2,5,i)), [], 'omitnan');
   
    %max theta, max beta, max alpha
    x(1,6,i) = max(thetas_temp, [], 'omitnan');
    x(3,6,i) = max(betas_temp(thetas_temp == x(1,6,i)), [], 'omitnan');
    x(2,6,i) = max(alphas_temp(thetas_temp == x(1,6,i) & betas_temp == x(3,6,i)), [], 'omitnan');
    
    %max theta, max beta, min alpha
    x(1,7,i) = max(thetas_temp, [], 'omitnan');
    x(3,7,i) = max(betas_temp(thetas_temp == x(1,7,i)), [], 'omitnan');
    x(2,7,i) = min(alphas_temp(thetas_temp == x(1,7,i) & betas_temp == x(3,7,i)), [], 'omitnan');
    
    %max theta, min beta, max alpha
    x(1,8,i) = max(thetas_temp, [], 'omitnan');
    x(3,8,i) = min(betas_temp(thetas_temp == x(1,8,i)), [], 'omitnan');
    x(2,8,i) = max(alphas_temp(thetas_temp == x(1,8,i) & betas_temp == x(3,8,i)), [], 'omitnan');
    
    %max theta, min beta, min alpha
    x(1,9,i) = max(thetas_temp, [], 'omitnan');
    x(3,9,i) = min(betas_temp(thetas_temp == x(1,9,i)), [], 'omitnan');
    x(2,9,i) = min(alphas_temp(thetas_temp == x(1,9,i) & betas_temp == x(3,9,i)), [], 'omitnan');
    
    % switch min theta
    
    %min theta, max alpha, max beta
    x(1,10,i) = min(thetas_temp, [], 'omitnan');
    x(2,10,i) = max(alphas_temp(thetas_temp == x(1,10,i)), [], 'omitnan');
    x(3,10,i) = max(betas_temp(thetas_temp == x(1,10,i) & alphas_temp == x(2,10,i)), [], 'omitnan');
    
    %min theta, max alpha, min beta
    x(1,11,i) = min(thetas_temp, [], 'omitnan');
    x(2,11,i) = max(alphas_temp(thetas_temp == x(1,1,i)), [], 'omitnan');
    el = min(betas_temp(thetas_temp == x(1,11,i) & alphas_temp == x(2,11,i)), [], 'omitnan');
    if ~isempty(el)
        x(3,11,i) = el;
    else
        x(3,11,i) = nan;
    end
    
    %min theta, min alpha, max beta
    x(1,12,i) = min(thetas_temp, [], 'omitnan');
    x(2,12,i) = min(alphas_temp(thetas_temp == x(1,12,i)), [], 'omitnan');
    el = max(betas_temp(thetas_temp == x(1,12,i) & alphas_temp == x(2,12,i)), [], 'omitnan');
    if ~isempty(el)
        x(3,12,i) = el;
    else
        x(3,12,i) = nan;
    end
        
    %min theta, min alpha, min beta
    x(1,13,i) = min(thetas_temp, [], 'omitnan');
    x(2,13,i) = min(alphas_temp(thetas_temp == x(1,13,i)), [], 'omitnan');
    el = min(betas_temp(thetas_temp == x(1,13,i) & alphas_temp == x(2,13,i)), [], 'omitnan');
    if ~isempty(el)
        x(3,13,i) = el;
    else
        x(3,13,i) = nan;
    end
   
    %min theta, max beta, max alpha
    x(1,14,i) = min(thetas_temp, [], 'omitnan');
    x(3,14,i) = max(betas_temp(thetas_temp == x(1,14,i)), [], 'omitnan');
    el = max(alphas_temp(thetas_temp == x(1,14,i) & betas_temp == x(3,14,i)), [], 'omitnan');
    if ~isempty(el)
        x(2,14,i) = el;
    else
        x(2,14,i) = nan;
    end
    
    %min theta, max beta, min alpha
    x(1,15,i) = min(thetas_temp, [], 'omitnan');
    x(3,15,i) = max(betas_temp(thetas_temp == x(1,15,i)), [], 'omitnan');
    el = min(alphas_temp(thetas_temp == x(1,15,i) & betas_temp == x(3,15,i)), [], 'omitnan');
    if ~isempty(el)
        x(2,15,i) = el;
    else
        x(2,15,i) = nan;
    end
    
    %min theta, min beta, max alpha
    x(1,16,i) = min(thetas_temp, [], 'omitnan');
    x(3,16,i) = min(betas_temp(thetas_temp == x(1,16,i)), [], 'omitnan');
    el = max(alphas_temp(thetas_temp == x(1,16,i) & betas_temp == x(3,16,i)), [], 'omitnan');
    if ~isempty(el)
        x(2,16,i) = el;
    else  
        x(2,16,i) =nan;
    end
    
    %min theta, min beta, min alpha
    x(1,17,i) = min(thetas_temp, [], 'omitnan');
    x(3,17,i) = min(betas_temp(thetas_temp == x(1,17,i)), [], 'omitnan');
    el = min(alphas_temp(thetas_temp == x(1,17,i) & betas_temp == x(3,17,i)), [], 'omitnan');
    if ~isempty(el)
       x(2,17,i) = el;
    else
        x(2,17,i) = nan;
    end
    
    % alpha lead %
    
    %max alpha, max theta, max beta
    x(2,18,i) = max(alphas_temp, [], 'omitnan');
    x(1,18,i) = max(thetas_temp(alphas_temp == x(2,18,i)), [], 'omitnan');
    el = max(betas_temp(alphas_temp == x(2,18,i) & thetas_temp == x(1,18,i)), [], 'omitnan');
    if ~isempty(el)
        x(3,18,i) = el;
    else
        x(3,18,i) = nan;
    end
    
    %max alpha, max theta, min beta
    x(2,19,i) = max(alphas_temp, [], 'omitnan');
    x(1,19,i) = max(thetas_temp(alphas_temp == x(2,19,i)), [], 'omitnan');
    el = min(betas_temp(alphas_temp == x(2,19,i) & thetas_temp == x(1,19,i)), [], 'omitnan');
    if ~isempty(el)
        x(3,19,i) = el;
    else
        x(3,19,i) = nan;
    end
    
    %max alpha, min theta, max beta
    x(2,20,i) = max(alphas_temp, [], 'omitnan');
    x(1,20,i) = min(thetas_temp(alphas_temp == x(2,20,i)), [], 'omitnan');
    el = max(betas_temp(alphas_temp == x(2,20,i) & thetas_temp == x(1,20,i)), [], 'omitnan');
    if ~isempty(el)
        x(3,20,i) = el;
    else
        x(3,20,i) = nan;
    end
  
    %max alpha, min theta, min beta
    x(2,21,i) = max(alphas_temp, [], 'omitnan');
    x(1,21,i) = min(thetas_temp(alphas_temp == x(2,21,i)), [], 'omitnan');
    el = min(betas_temp(alphas_temp == x(2,21,i) & thetas_temp == x(1,21,i)), [], 'omitnan');
    if ~isempty(el)
        x(3,21,i) = el;
    else
        x(3,21,i) = nan;
    end
    
    %max alpha, max beta, max theta
    x(2,22,i) = max(alphas_temp, [], 'omitnan');
    x(3,22,i) = max(betas_temp(alphas_temp == x(2,22,i)), [], 'omitnan');
    el = max(thetas_temp(alphas_temp == x(2,22,i) & betas_temp == x(3,22,i)), [], 'omitnan');
    if ~isempty(el)
        x(1,22,i) = el;
    else
        x(1,22,i) = nan;
    end
    
    %max alpha, max beta, min theta
    x(2,23,i) = max(alphas_temp, [], 'omitnan');
    x(3,23,i) = max(betas_temp(alphas_temp == x(2,23,i)), [], 'omitnan');
    el = min(thetas_temp(alphas_temp == x(2,23,i) & betas_temp == x(3,23,i)), [], 'omitnan');
    if ~isempty(el)
        x(1,23,i) = el;
    else
        x(1,23,i) = nan;
    end
    
    %max alpha, min beta, max theta
    x(2,24,i) = max(alphas_temp, [], 'omitnan');
    x(3,24,i) = min(betas_temp(alphas_temp == x(2,24,i)), [], 'omitnan');
    el = max(thetas_temp(alphas_temp == x(2,24,i) & betas_temp == x(3,24,i)), [], 'omitnan');
    if ~isempty(el)
        x(1,24,i) = el;
    else
        x(1,24,i) = nan;
    end
    
    %max alpha, min beta, min theta
    x(2,25,i) = max(alphas_temp, [], 'omitnan');
    x(3,25,i) = min(betas_temp(alphas_temp == x(2,25,i)), [], 'omitnan');
    el =  min(thetas_temp(alphas_temp == x(2,25,i) & betas_temp == x(3,25,i)), [], 'omitnan');
    if ~isempty(el)
        x(1,25,i) = el;
    else
        x(1,25,i) = nan;
    end
    
    % switch min alpha
    
    %min alpha, max theta, max beta
    x(2,26,i) = min(alphas_temp, [], 'omitnan');
    x(1,26,i) = max(thetas_temp(alphas_temp == x(2,26,i)), [], 'omitnan');
    x(3,26,i) = max(betas_temp(alphas_temp == x(2,26,i) & thetas_temp == x(1,26,i)), [], 'omitnan');
    
    %min alpha, max theta, min beta
    x(2,27,i) = min(alphas_temp, [], 'omitnan');
    x(1,27,i) = max(thetas_temp(alphas_temp == x(2,27,i)), [], 'omitnan');
    x(3,27,i) = min(betas_temp(alphas_temp == x(2,27,i) & thetas_temp == x(1,27,i)), [], 'omitnan');
    
    %min alpha, min theta, max beta
    x(2,28,i) = min(alphas_temp, [], 'omitnan');
    x(1,28,i) = min(thetas_temp(alphas_temp == x(2,28,i)), [], 'omitnan');
    x(3,28,i) = max(betas_temp(alphas_temp == x(2,28,i) & thetas_temp == x(1,28,i)), [], 'omitnan');
    
    %min alpha, min theta, min beta
    x(2,29,i) = min(alphas_temp, [], 'omitnan');
    x(1,29,i) = min(thetas_temp(alphas_temp == x(2,29,i)), [], 'omitnan');
    x(3,29,i) = min(betas_temp(alphas_temp == x(2,29,i) & thetas_temp == x(1,29,i)), [], 'omitnan');
    
    %min alpha, max beta, max theta
    x(2,30,i) = min(alphas_temp, [], 'omitnan');
    x(3,30,i) = max(betas_temp(alphas_temp == x(2,30,i)), [], 'omitnan');
    x(1,30,i) = max(thetas_temp(alphas_temp == x(2,30,i) & betas_temp == x(3,30,i)), [], 'omitnan');
    
    %min alpha, max beta, min theta
    x(2,31,i) = min(alphas_temp, [], 'omitnan');
    x(3,31,i) = max(betas_temp(alphas_temp == x(2,31,i)), [], 'omitnan');
    x(1,31,i) = min(thetas_temp(alphas_temp == x(2,31,i) & betas_temp == x(3,31,i)), [], 'omitnan');
    
    %min alpha, min beta, max theta
    x(2,32,i) = min(alphas_temp, [], 'omitnan');
    x(3,32,i) = min(betas_temp(alphas_temp == x(2,32,i)), [], 'omitnan');
    el = max(thetas_temp(alphas_temp == x(2,32,i) & betas_temp == x(3,33,i)), [], 'omitnan');
    if ~isnan(el)
        x(1,32,i) = el;
    else
        x(1,32,i) = nan;
    end
    
    %min alpha, min beta, min theta
    x(2,33,i) = min(alphas_temp, [], 'omitnan');
    x(3,33,i) = min(betas_temp(alphas_temp == x(2,33,i)), [], 'omitnan');
    x(1,33,i) = min(thetas_temp(alphas_temp == x(2,33,i) & betas_temp == x(3,33,i)), [], 'omitnan');
    
    
    % beta lead %
    
   
    %max beta, max theta, max alpha
    x(3,34,i) = max(betas_temp, [], 'omitnan');
    x(1,34,i) = max(thetas_temp(betas_temp == x(3,34,i)), [], 'omitnan');
    x(2,34,i) = max(alphas_temp(betas_temp == x(3,34,i) & thetas_temp == x(1,34,i)), [], 'omitnan');
    
    %max beta, max theta, min alpha
    x(3,35,i) = max(betas_temp, [], 'omitnan');
    x(1,35,i) = max(thetas_temp(betas_temp == x(3,35,i)), [], 'omitnan');
    x(2,35,i) = min(alphas_temp(betas_temp == x(3,35,i) & thetas_temp == x(1,35,i)), [], 'omitnan');
    
    %max beta, min theta, max alpha
    x(3,36,i) = max(betas_temp, [], 'omitnan');
    x(1,36,i) = min(thetas_temp(betas_temp == x(3,36,i)), [], 'omitnan');
    x(2,36,i) = max(alphas_temp(betas_temp == x(3,36,i) & thetas_temp == x(1,36,i)), [], 'omitnan');
    
    %max beta, min theta, min alpha
    x(3,37,i) = max(betas_temp, [], 'omitnan');
    x(1,37,i) = min(thetas_temp(betas_temp == x(3,37,i)), [], 'omitnan');
    x(2,37,i) = min(alphas_temp(betas_temp == x(3,37,i) & thetas_temp == x(1,37,i)), [], 'omitnan');
    
    %max beta, max alpha, max theta
    x(3,38,i) = max(betas_temp, [], 'omitnan');
    x(2,38,i) = max(alphas_temp(betas_temp == x(3,38,i)), [], 'omitnan');
    x(1,38,i) = max(thetas_temp(betas_temp == x(3,38,i) & alphas_temp == x(2,38,i)), [], 'omitnan');
    
    %max beta, max alpha, min theta
    x(3,39,i) = max(betas_temp, [], 'omitnan');
    x(2,39,i) = max(alphas_temp(betas_temp == x(3,39,i)), [], 'omitnan');
    x(1,39,i) = min(thetas_temp(betas_temp == x(3,39,i) & alphas_temp == x(2,39,i)), [], 'omitnan');
    
    %max beta, min alpha, max theta
    x(3,40,i) = max(betas_temp, [], 'omitnan');
    x(2,40,i) = min(alphas_temp(betas_temp == x(3,40,i)), [], 'omitnan');
    x(1,40,i) = max(thetas_temp(betas_temp == x(3,40,i) & alphas_temp == x(2,40,i)), [], 'omitnan');
    
    %max beta, min alpha, min theta
    x(3,41,i) = max(betas_temp, [], 'omitnan');
    x(2,41,i) = min(alphas_temp(betas_temp == x(3,41,i)), [], 'omitnan');
    x(1,41,i) = min(thetas_temp(betas_temp == x(3,41,i) & alphas_temp == x(2,41,i)), [], 'omitnan');
    
   % switch min beta
    
   %min beta, max theta, max alpha
    x(3,42,i) = min(betas_temp, [], 'omitnan');
    x(1,42,i) = max(thetas_temp(betas_temp == x(3,42,i)), [], 'omitnan');
    x(2,42,i) = max(alphas_temp(betas_temp == x(3,42,i) & thetas_temp == x(1,42,i)), [], 'omitnan');
    
    %min beta, max theta, min alpha
    x(3,43,i) = min(betas_temp, [], 'omitnan');
    x(1,43,i) = max(thetas_temp(betas_temp == x(3,43,i)), [], 'omitnan');
    x(2,43,i) = min(alphas_temp(betas_temp == x(3,43,i) & thetas_temp == x(1,43,i)), [], 'omitnan');
    
    %min beta, min theta, max alpha
    x(3,44,i) = min(betas_temp, [], 'omitnan');
    x(1,44,i) = min(thetas_temp(betas_temp == x(3,44,i)), [], 'omitnan');
    x(2,44,i) = max(alphas_temp(betas_temp == x(3,44,i) & thetas_temp == x(1,44,i)), [], 'omitnan');
    
    %min beta, min theta, min alpha
    x(3,45,i) = min(betas_temp, [], 'omitnan');
    x(1,45,i) = min(thetas_temp(betas_temp == x(3,45,i)), [], 'omitnan');
    x(2,45,i) = min(alphas_temp(betas_temp == x(3,45,i) & thetas_temp == x(1,45,i)), [], 'omitnan');
    
    %min beta, max alpha, max theta
    x(3,46,i) = min(betas_temp, [], 'omitnan');
    x(2,46,i) = max(alphas_temp(betas_temp == x(3,46,i)), [], 'omitnan');
    x(1,46,i) = max(thetas_temp(betas_temp == x(3,46,i) & alphas_temp == x(2,46,i)), [], 'omitnan');
    
    %min beta, max alpha, min theta
    x(3,47,i) = min(betas_temp, [], 'omitnan');
    x(2,47,i) = max(alphas_temp(betas_temp == x(3,47,i)), [], 'omitnan');
    x(1,47,i) = min(thetas_temp(betas_temp == x(3,47,i) & alphas_temp == x(2,47,i)), [], 'omitnan');
    
    %min beta, min alpha, max theta
    x(3,48,i) = min(betas_temp, [], 'omitnan');
    x(2,48,i) = min(alphas_temp(betas_temp == x(3,48,i)), [], 'omitnan');
    x(1,48,i) = max(thetas_temp(betas_temp == x(3,48,i) & alphas_temp == x(2,48,i)), [], 'omitnan');
    
    %min beta, min alpha, min theta
    x(3,49,i) = min(betas_temp, [], 'omitnan');
    x(2,49,i) = min(alphas_temp(betas_temp == x(3,49,i)), [], 'omitnan');
    x(1,49,i) = min(thetas_temp(betas_temp == x(3,49,i) & alphas_temp == x(2,49,i)), [], 'omitnan');
   
    i = i + 1;
end
 
nsamples = 10000;
 
%% Analyis
 
% generate distribution data
[dist.z,dist.target] = dist.generate(dist);

itrial = randi(dist.ntrials,[nsamples 1]);
 %power
scale = @(p, alpha, beta) 1 + (1/nchoosek(size(p, 2), 3)) .* alpha .* sum( (max((reshape(p(:, nchoosek(1:size(p, 2), 3)), [size(p, 1), nchoosek(size(p, 2), 3), 3])), [], 3) - min((reshape(p(:, nchoosek(1:size(p, 2), 3)), [size(p, 1), nchoosek(size(p, 2), 3), 3])), [], 3)).^(beta) .* prod( reshape(p(:, nchoosek(1:size(p, 2), 3)), [size(p, 1), nchoosek(size(p, 2), 3), 3]), 3), 2);
 
% normalize weights
c = max(wcs, [], 2); ws = wcs./c;
 
% threshold crossing
[e, rt] = deal(nan(nsamples,1));
for isample = 1:nsamples
    for combos = 1:size(x, 2)
        for wces = 1:size(x,3)
            sc = scale( exp( dist.z( 2:end, :, itrial(isample) ) ), x(2,combos,wces), x(3,combos,wces) );
            sc = [sc(1:end-1); 0]; % solve crossing problem
            isdec = exp(dist.z(2:end, :, itrial(isample))) > x(1,combos,wces)*sc.*ones(1,N);
            if ~isempty(find( sum( isdec ,2), 1))
                rt(isample,combos,wces) = find( sum( isdec ,2), 1);
                [~, n] = max( dist.z(rt(isample,combos,wces)+1, :, itrial(isample)) );
                sc = scale( exp( dist.z( find( sum( isdec ,2), 1), :, itrial(isample) ) ), x(2,combos,wces), x(3,combos,wces) );
                bound(isample,combos,wces) = x(1,combos,wces)*sc;
                e(isample,combos,wces) = n~=dist.target(itrial(isample));
            else
                rt(isample,combos,wces) = nan;
                bound(isample,combos,wces) = nan;
                e(isample,combos,wces) = nan;
            end
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
            bound_avg{combos,wces}(r_t)=mean(bound(rnk==r_t,combos,wces));
            bound_var{combos,wces}(r_t)=var(bound(rnk==r_t,combos,wces));
            bound_max{combos,wces}(r_t)=max(bound(rnk==r_t,combos,wces));
            bound_min{combos,wces}(r_t)=min(bound(rnk==r_t,combos,wces));
            
            log_bound = log(bound(rnk==r_t,combos,wces)./(1-bound(rnk==r_t,combos,wces)));
            
            bound_avg_logd{combos,wces}(r_t)=mean(log_bound);
            bound_var_logd{combos,wces}(r_t)=var(log_bound);
            bound_max_logd{combos,wces}(r_t)=max(log_bound);
            bound_min_logd{combos,wces}(r_t)=min(log_bound);
            
       end
   end
end

%% power

N = 3;
C = [30, 136, 229; 216, 27, 96; 255, 193, 7]./255;

i=1;
%flat
combos_i = [50, 58, 50, 52, 36, 63, 64];
wces = [1, 1, 2, 2, 3, 3, 3];
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
hf = gcf;
ha = gca;
set(hf, 'color', 'w', 'units', 'inches', 'position', [0 0 4 4])
set(ha, 'color', 'w', 'units', 'inches', 'position', [1 1 2.75 2.5])
ha.XAxis.TickValues = [0 10 20 30];
ha.YAxis.TickValues = [0 2 4 6];
xlabel('reaction time ({\it rT})')

if ~exist('figs', 'dir')
       mkdir('figs')
    end

export_fig(append('./figs/fig8h'), '-pdf', '-eps', '-q101');
savefig([pwd '/figs/fig8h'])

%%

i=1;
%decreasing
combos_i = [53, 54, 53];
wces = [1, 1, 2];
figure
for combos = combos_i
    scatter(rt_b_avg{combos,wces(i)}(2:end), bound_avg_logd{combos,wces(i)}(2:end), 10, C(wces(i),:), 'filled')
    hold on
    
    gprMdl_mn = fitrgp(rt_b_avg{combos,wces(i)}(2:end), real(bound_avg_logd{combos,wces(i)}(2:end)),'Basis','pureQuadratic','FitMethod','exact',...
        'PredictMethod','exact','SigmaLowerBound',0.251);
    ypred = resubPredict(gprMdl_mn);
    X = rt_b_avg{combos,wces(i)}(~isnan(real(bound_avg_logd{combos,wces(i)})));
    if size(X)==size(ypred)
        plot(X, ypred, 'Color', C(wces(i),:), 'LineWidth', 2)
    else
        plot(X(2:end), ypred, 'Color', C(wces(i),:), 'LineWidth', 2)
    end
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

if ~exist('figs', 'dir')
       mkdir('figs')
    end

export_fig(append('./figs/fig8e'), '-pdf', '-eps', '-q101');
savefig([pwd '/figs/fig8e'])

%%

i=1;
%increasing
combos_i = [3, 21, 2, 3, 4, 3, 52, 24];
wces = [1, 1, 2, 2, 3, 3, 3, 3];
figure
for combos = combos_i
    scatter(rt_b_avg{combos,wces(i)}(2:end), bound_avg_logd{combos,wces(i)}(2:end), 10, C(wces(i),:), 'filled')
    hold on
    
    gprMdl_mn = fitrgp(rt_b_avg{combos,wces(i)}(2:end), real(bound_avg_logd{combos,wces(i)}(2:end)),'Basis','pureQuadratic','FitMethod','exact',...
        'PredictMethod','exact','SigmaLowerBound',0.1);
    ypred = resubPredict(gprMdl_mn);
    X = rt_b_avg{combos,wces(i)}(~isnan(real(bound_avg_logd{combos,wces(i)})));
    if size(X)==size(ypred)
        plot(X, ypred, 'Color', C(wces(i),:), 'LineWidth', 2)
    else
        plot(X(2:end), ypred, 'Color', C(wces(i),:), 'LineWidth', 2)
    end
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
title(ha, 'power(\theta, \alpha, \beta)', 'FontSize', 20);
ha.XAxis.TickValues = [0 10 20 30];
ha.YAxis.TickValues = [0 2 4 6];

if ~exist('figs', 'dir')
       mkdir('figs')
    end

export_fig(append('./figs/fig8b'), '-pdf', '-eps', '-q101');
savefig([pwd '/figs/fig8b'])