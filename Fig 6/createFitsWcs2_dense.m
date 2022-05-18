function [datx, daty, fitresult, fts, optss, thetas, alphas, r_opt] = createFitsWcs2_dense(figs, figs2)
%CREATEFITS(X_AXIS_1,Y_AXIS_1,X_AXIS_2,Y_AXIS_2,X_AXIS_3,Y_AXIS_3,X_AXIS_4,Y_AXIS_4,X_AXIS_5,Y_AXIS_5)
%  Create fits.
%
%  Data for 'Wcs2_alph1' fit:
%      X Input : x_axis_1
%      Y Output: y_axis_1
%  Data for 'Wcs2_alph2' fit:
%      X Input : x_axis_2
%      Y Output: y_axis_2
%  Data for 'Wcs2_alph3' fit:
%      X Input : x_axis_3
%      Y Output: y_axis_3
%  Data for 'Wcs2_alph4' fit:
%      X Input : x_axis_4
%      Y Output: y_axis_4
%  Data for 'Wcs2_alph5' fit:
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

wcsss = 2;
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

%% Fit: 'Wcs2_alph1'.
[xData, yData] = prepareCurveData( x_axis_1, y_axis_1 );

% Set up fittype and options.
ft = fittype( 'gauss2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 -Inf -Inf 0];
opts.StartPoint = [-0.383807656565691 0.658470637957963 0.0774282937910058 -0.383807656565691 0.735898931748969 0.0774282937910058];

% Fit model to data.
[fitresult{1}, gof(1)] = fit( xData, yData, ft, opts );

if figs
    % Plot fit with data.
    figure( 'Name', 'Wcs2_alph1' );
    h = plot( fitresult{1}, xData, yData, 'predobs' );
    legend( h, 'y_axis_1 vs. x_axis_1', 'Wcs2_alph1', 'Lower bounds (Wcs2_alph1)', 'Upper bounds (Wcs2_alph1)', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % Label axes
    xlabel( 'x_axis_1', 'Interpreter', 'none' );
    ylabel( 'y_axis_1', 'Interpreter', 'none' );
    grid on
end

datx{1}=xData;
daty{1}=yData;
fts{1}=ft;
optss{1}=opts;

%% Fit: 'Wcs2_alph2'.
[xData, yData] = prepareCurveData( x_axis_2, y_axis_2 );

% Set up fittype and options.
ft = fittype( 'gauss2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 -Inf -Inf 0];
opts.StartPoint = [-0.383631636363678 0.681553007150304 0.0660659262658488 -0.383631636363678 0.747618933416153 0.0660659262658488];

% Fit model to data.
[fitresult{2}, gof(2)] = fit( xData, yData, ft, opts );

if figs
    % Plot fit with data.
    figure( 'Name', 'Wcs2_alph2' );
    h = plot( fitresult{2}, xData, yData, 'predobs' );
    legend( h, 'y_axis_2 vs. x_axis_2', 'Wcs2_alph2', 'Lower bounds (Wcs2_alph2)', 'Upper bounds (Wcs2_alph2)', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % Label axes
    xlabel( 'x_axis_2', 'Interpreter', 'none' );
    ylabel( 'y_axis_2', 'Interpreter', 'none' );
    grid on
end

datx{2}=xData;
daty{2}=yData;
fts{2}=ft;
optss{2}=opts;

%% Fit: 'Wcs2_alph3'.
[xData, yData] = prepareCurveData( x_axis_3, y_axis_3 );

% Set up fittype and options.
ft = fittype( 'gauss2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 -Inf -Inf 0];
opts.StartPoint = [-0.383544282828304 0.712351861774209 0.0546204290143528 -0.383544282828304 0.766972290788562 0.0546204290143528];

% Fit model to data.
[fitresult{3}, gof(3)] = fit( xData, yData, ft, opts );

if figs
    % Plot fit with data.
    figure( 'Name', 'Wcs2_alph3' );
    h = plot( fitresult{3}, xData, yData, 'predobs' );
    legend( h, 'y_axis_3 vs. x_axis_3', 'Wcs2_alph3', 'Lower bounds (Wcs2_alph3)', 'Upper bounds (Wcs2_alph3)', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % Label axes
    xlabel( 'x_axis_3', 'Interpreter', 'none' );
    ylabel( 'y_axis_3', 'Interpreter', 'none' );
    grid on
end

datx{3}=xData;
daty{3}=yData;
fts{3}=ft;
optss{3}=opts;

%% Fit: 'Wcs2_alph4'.
[xData, yData] = prepareCurveData( x_axis_4, y_axis_4 );

% Set up fittype and options.
ft = fittype( 'gauss2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 -Inf -Inf 0];
opts.StartPoint = [-0.382854000000029 0.738405871079533 0.0441295454091938 -0.382854000000029 0.782535416488727 0.0441295454091938];

% Fit model to data.
[fitresult{4}, gof(4)] = fit( xData, yData, ft, opts );

if figs
    % Plot fit with data.
    figure( 'Name', 'Wcs2_alph4' );
    h = plot( fitresult{4}, xData, yData, 'predobs' );
    legend( h, 'y_axis_4 vs. x_axis_4', 'Wcs2_alph4', 'Lower bounds (Wcs2_alph4)', 'Upper bounds (Wcs2_alph4)', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % Label axes
    xlabel( 'x_axis_4', 'Interpreter', 'none' );
    ylabel( 'y_axis_4', 'Interpreter', 'none' );
    grid on
end

datx{4}=xData;
daty{4}=yData;
fts{4}=ft;
optss{4}=opts;

%% Fit: 'Wcs2_alph5'.
[xData, yData] = prepareCurveData( x_axis_5, y_axis_5 );

% Set up fittype and options.
ft = fittype( 'gauss2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 -Inf -Inf 0];
opts.StartPoint = [-0.382553797979813 0.762351237207253 0.0378284169989752 -0.382553797979813 0.800179654206229 0.0378284169989752];

% Fit model to data.
[fitresult{5}, gof(5)] = fit( xData, yData, ft, opts );

if figs
    % Plot fit with data.
    figure( 'Name', 'Wcs2_alph5' );
    h = plot( fitresult{5}, xData, yData, 'predobs' );
    legend( h, 'y_axis_5 vs. x_axis_5', 'Wcs2_alph5', 'Lower bounds (Wcs2_alph5)', 'Upper bounds (Wcs2_alph5)', 'Location', 'NorthEast', 'Interpreter', 'none' );
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
        plot( fitresult{i}, datx{i}, daty{i}, 'predobs' );
    end
end
