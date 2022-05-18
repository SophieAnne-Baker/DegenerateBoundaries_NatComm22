function [datx, daty, fitresult, fts, optss, thetas, alphas, r_opt] = createFitsWcs1_dense(figs, figs2)
%CREATEFITS(X_AXIS_1,Y_AXIS_1,X_AXIS_2,Y_AXIS_2,X_AXIS_3,Y_AXIS_3,X_AXIS_4,Y_AXIS_4,X_AXIS_5,Y_AXIS_5)
%  Create fits.
%
%  Data for 'Wcs1_alph1' fit:
%      X Input : x_axis_1
%      Y Output: y_axis_1
%  Data for 'Wcs1_alph2' fit:
%      X Input : x_axis_2
%      Y Output: y_axis_2
%  Data for 'Wcs1_alph3' fit:
%      X Input : x_axis_3
%      Y Output: y_axis_3
%  Data for 'Wcs1_alph4' fit:
%      X Input : x_axis_4
%      Y Output: y_axis_4
%  Data for 'Wcs1_alph5' fit:
%      X Input : x_axis_5
%      Y Output: y_axis_5
%  Output:
%      fitresult : a cell-array of fit objects representing the fits.
%      gof : structure array with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.


%% data

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

wcsss = 1;
x_axis_1 = thetas(1,:,wcsss);
x_axis_2 = thetas(2,:,wcsss);
x_axis_3 = thetas(3,:,wcsss);
x_axis_4 = thetas(4,:,wcsss);
x_axis_5 = thetas(5,:,wcsss);
y_axis_1 = r_opt(1,:,wcsss);
y_axis_2 = r_opt(2,:,wcsss);
y_axis_3 = r_opt(3,:,wcsss);
y_axis_4 = r_opt(4,:,wcsss);
y_axis_5 = r_opt(5,:,wcsss);


%% Initialization.

% Initialize arrays to store fits and goodness-of-fit.
fitresult = cell( 5, 1 );
gof = struct( 'sse', cell( 5, 1 ), ...
    'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );

%% Fit: 'Wcs1_alph1'.
[xData, yData] = prepareCurveData( x_axis_1, y_axis_1 );

% Set up fittype and options.
ft = fittype( 'gauss2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 -Inf -Inf 0];
opts.StartPoint = [-0.130505313131256 0.93335147506079 0.0207364246805322 -0.130505313131256 0.954087899741322 0.0207364246805322];

% Fit model to data.
[fitresult{1}, gof(1)] = fit( xData, yData, ft, opts );

if figs
    % Plot fit with data.
    figure( 'Name', 'Wcs1_alph1' );
    h = plot( fitresult{1}, xData, yData, 'predobs' );
    legend( h, 'y_axis_1 vs. x_axis_1', 'Wcs1_alph1', 'Lower bounds (Wcs1_alph1)', 'Upper bounds (Wcs1_alph1)', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % Label axes
    xlabel( 'x_axis_1', 'Interpreter', 'none' );
    ylabel( 'y_axis_1', 'Interpreter', 'none' );
    grid on
end

datx{1}=xData;
daty{1}=yData;
fts{1}=ft;
optss{1}=opts;

%% Fit: 'Wcs1_alph2'.
[xData, yData] = prepareCurveData( x_axis_2, y_axis_2 );

% Set up fittype and options.
ft = fittype( 'gauss2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 -Inf -Inf 0];
opts.StartPoint = [-0.130598484848428 0.93666330718391 0.0180103027265913 -0.130598484848428 0.954673609910501 0.0180103027265913];

% Fit model to data.
[fitresult{2}, gof(2)] = fit( xData, yData, ft, opts );

if figs
    % Plot fit with data.
    figure( 'Name', 'Wcs1_alph2' );
    h = plot( fitresult{2}, xData, yData, 'predobs' );
    legend( h, 'y_axis_2 vs. x_axis_2', 'Wcs1_alph2', 'Lower bounds (Wcs1_alph2)', 'Upper bounds (Wcs1_alph2)', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % Label axes
    xlabel( 'x_axis_2', 'Interpreter', 'none' );
    ylabel( 'y_axis_2', 'Interpreter', 'none' );
    grid on
end

datx{2}=xData;
daty{2}=yData;
fts{2}=ft;
optss{2}=opts;

%% Fit: 'Wcs1_alph3'.
[xData, yData] = prepareCurveData( x_axis_3, y_axis_3 );

% Set up fittype and options.
ft = fittype( 'gauss2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 -Inf -Inf 0];
opts.StartPoint = [-0.129902101010047 0.938978779240275 0.0174944079223418 -0.129902101010047 0.956473187162617 0.0174944079223418];

% Fit model to data.
[fitresult{3}, gof(3)] = fit( xData, yData, ft, opts );

if figs
    % Plot fit with data.
    figure( 'Name', 'Wcs1_alph3' );
    h = plot( fitresult{3}, xData, yData, 'predobs' );
    legend( h, 'y_axis_3 vs. x_axis_3', 'Wcs1_alph3', 'Lower bounds (Wcs1_alph3)', 'Upper bounds (Wcs1_alph3)', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % Label axes
    xlabel( 'x_axis_3', 'Interpreter', 'none' );
    ylabel( 'y_axis_3', 'Interpreter', 'none' );
    grid on
end

datx{3}=xData;
daty{3}=yData;
fts{3}=ft;
optss{3}=opts;

%% Fit: 'Wcs1_alph4'.
[xData, yData] = prepareCurveData( x_axis_4, y_axis_4 );

% Set up fittype and options.
ft = fittype( 'gauss2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 -Inf -Inf 0];
opts.StartPoint = [-0.130375353535299 0.937629916709574 0.0181688391876924 -0.130375353535299 0.955798755897267 0.0181688391876924];

% Fit model to data.
[fitresult{4}, gof(4)] = fit( xData, yData, ft, opts );

if figs
    % Plot fit with data.
    figure( 'Name', 'Wcs1_alph4' );
    h = plot( fitresult{4}, xData, yData, 'predobs' );
    legend( h, 'y_axis_4 vs. x_axis_4', 'Wcs1_alph4', 'Lower bounds (Wcs1_alph4)', 'Upper bounds (Wcs1_alph4)', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % Label axes
    xlabel( 'x_axis_4', 'Interpreter', 'none' );
    ylabel( 'y_axis_4', 'Interpreter', 'none' );
    grid on
end

datx{4}=xData;
daty{4}=yData;
fts{4}=ft;
optss{4}=opts;

%% Fit: 'Wcs1_alph5'.
[xData, yData] = prepareCurveData( x_axis_5, y_axis_5 );

% Set up fittype and options.
ft = fittype( 'gauss2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 -Inf -Inf 0];
opts.StartPoint = [-0.130329414141356 0.938708769486822 0.0176294127990685 -0.130329414141356 0.956338182285891 0.0176294127990685];

% Fit model to data.
[fitresult{5}, gof(5)] = fit( xData, yData, ft, opts );

if figs
    % Plot fit with data.
    figure( 'Name', 'Wcs1_alph5' );
    h = plot( fitresult{5}, xData, yData, 'predobs' );
    legend( h, 'y_axis_5 vs. x_axis_5', 'Wcs1_alph5', 'Lower bounds (Wcs1_alph5)', 'Upper bounds (Wcs1_alph5)', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % Label axes
    xlabel( 'x_axis_5', 'Interpreter', 'none' );
    ylabel( 'y_axis_5', 'Interpreter', 'none' );
    grid on
end

datx{5}=xData;
daty{5}=yData;
fts{5}=ft;
optss{5}=opts;

%% combined

if figs2
    figure; hold on;
    for i=1:5
        plot( fitresult{i}, thetas(i,:,wcsss), r_opt(i,:,wcsss), 'predobs' );
    end
end

