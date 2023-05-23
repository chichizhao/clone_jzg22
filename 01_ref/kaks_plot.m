%%/bin/matlab
% author: zhaochichi github: https://github.com/chichizhao/clone_jzg22
% this script is count the ks distrubution number of the ks 
% 5 species ananas oryza sorghum Ta and Tl
% first we need get each ks value from test
ananasKs = test.Ksananas;
oryzaKs = test.Ksoryza;
sorghumKs = test.Kssorghum;
TaKs = test.KsTa;

% we need remove some abnormal value
ananasKs = ananasKs(ananasKs<8);
oryzaKs = oryzaKs(oryzaKs<8);
sorghumKs = sorghumKs(sorghumKs<8);
TaKs = TaKs(TaKs<8);

% we set the distance of each group is 0.1
distance = 0.1;
% we conut distribution of each group ks from 0 to 8 in 0.1 distance
ananasKsDistribution = hist(ananasKs,0:distance:8);
oryzaKsDistribution = hist(oryzaKs,0:distance:8);
sorghumKsDistribution = hist(sorghumKs,0:distance:8);
TaKsDistribution = hist(TaKs,0:distance:8);
% we need to normalize the distribution
ananasKsDistribution = ananasKsDistribution/sum(ananasKsDistribution);
oryzaKsDistribution = oryzaKsDistribution/sum(oryzaKsDistribution);
sorghumKsDistribution = sorghumKsDistribution/sum(sorghumKsDistribution);
TaKsDistribution = TaKsDistribution/sum(TaKsDistribution);
% we need to plot the distribution in a plot that contain 5 sub plots
x=0:0.1:8
subplot(5,1,1);
[xData, yData] = prepareCurveData( x, ananasKsDistribution );
% Set up fittype and options.
ft = fittype( 'gauss2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 -Inf -Inf 0];
opts.Normalize = 'on';
opts.Robust = 'LAR';
opts.StartPoint = [0.0587920897915553 -1.06262785510556 0.131880802117382 0.0392228358843854 -1.3176585403309 0.252169313757329];
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
h = plot( fitresult, xData, yData, 'predobs', 0.9 );
legend( h, 'Ananas', ' Gaussian', ' - 10 %', ' + 10%', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'Ks', 'Interpreter', 'none' );
ylabel( 'Density', 'Interpreter', 'none' );
grid on


subplot(5,1,2);
[xData, yData] = prepareCurveData( x, oryzaKsDistribution );

% Set up fittype and options.
ft = fittype( 'gauss3' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0];
opts.Normalize = 'on';
opts.Robust = 'LAR';
opts.StartPoint = [0.0609037328094303 -1.36016365453512 0.112627994339289 0.0381455760746419 -1.10513296930979 0.146005537046995 0.0349705304518664 -0.17002045681689 0.124503878474398];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
h = plot( fitresult, xData, yData, 'predobs', 0.9 );
legend( h, 'Oryza', ' Gaussian', ' - 10 %', ' + 10%', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'Ks', 'Interpreter', 'none' );
ylabel( 'Density', 'Interpreter', 'none' );
grid on


subplot(5,1,3)
[xData, yData] = prepareCurveData( x, sorghumKsDistribution );

% Set up fittype and options.
ft = fittype( 'gauss3' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0];
opts.Normalize = 'on';
opts.Robust = 'LAR';
opts.StartPoint = [0.0489370236662655 -1.3176585403309 0.128179134970366 0.0373044524669073 0 0.143454851405501 0.0300412599673342 -0.977617626697119 0.148180084086096];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
h = plot( fitresult, xData, yData, 'predobs', 0.9 );
legend( h, ' Sorghum', ' Gaussian', ' - 10 %', ' + 10%', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'Ks', 'Interpreter', 'none' );
ylabel( 'Density', 'Interpreter', 'none' );
grid on
subplot(5,1,4)
[xData, yData] = prepareCurveData( x, TaKsDistribution );

% Set up fittype and options.
ft = fittype( 'gauss3' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0];
opts.Normalize = 'on';
opts.Robust = 'LAR';
opts.StartPoint = [0.0702357563850688 -1.36016365453512 0.103103313107128 0.0525370522498002 -1.06262785510556 0.101621846128498 0.0260265120572487 -0.892607398288673 0.151830629567967];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
h = plot( fitresult, xData, yData, 'predobs', 0.9 );
legend( h, ' Typha angustifolia', ' Gaussian', ' - 10 %', ' + 10%', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'Ks', 'Interpreter', 'none' );
ylabel( 'Density', 'Interpreter', 'none' );
grid on