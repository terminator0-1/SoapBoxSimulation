clear; clc; close all;

%% Figure 1: Kart Speed Verification
T1 = readtable('KartSpeedVerificaton.csv', 'TextType', 'string');

% Separate the series (cut simulated curves at 23 m; keep GPS full so it reaches the 23 m edge)
dyn1  = T1(T1.series == "Dynamic"             & T1.x_m <= 23, :);
eb1   = T1(T1.series == "Energy balance"      & T1.x_m <= 23, :);
gps1  = T1(T1.series == "GPS Speed"           & T1.x_m <= 23, :);
dBand = T1(T1.series == "Dynamic +/-15% band" & T1.x_m <= 23, :);
gBand = T1(T1.series == "GPS +/-15% band"     & T1.x_m <= 23, :);

figure(1);
hold on; grid on;

% Plot lines (change colors/styles right inside plot: 'b', 'r', 'g', 'k', etc.)
plot(dyn1.x_m, dyn1.speed_mps, 'b-',  'LineWidth', 1.5);
plot(eb1.x_m,  eb1.speed_mps,  'g--', 'LineWidth', 1.5);
plot(gps1.x_m, gps1.speed_mps, 'k',   'LineWidth', 1.5);

% +/- 15% bounds
% plot(dBand.x_m, dBand.lower_speed_mps, 'b:');
% plot(dBand.x_m, dBand.upper_speed_mps, 'b:', 'HandleVisibility', 'off');
% 
% plot(gBand.x_m, gBand.lower_speed_mps, 'g:');
% plot(gBand.x_m, gBand.upper_speed_mps, 'g:', 'HandleVisibility', 'off');

xlim([0 23]);
xlabel('Distance (m)');
ylabel('Speed (m/s)');
title('Kart Speed Verification');
legend('Dynamic', 'Energy balance', 'GPS Speed', 'Location', 'best');
hold off;


%% Figure 2: Course Speed Comparison
T2 = readtable('CourseSpeedComparison.csv', 'TextType', 'string');

% Cut simulated at 109 m; include the 113.93 m point for Real kart so the line reaches 113 m
dyn2   = T2(T2.series == "Dynamic"               & T2.x_m <= 109, :);
eb2    = T2(T2.series == "Energy balance"        & T2.x_m <= 109, :);
real2  = T2(T2.series == "Real kart speed"       & T2.x_m <= 114, :);
dBand2 = T2(T2.series == "Dynamic +/-15% band"   & T2.x_m <= 109, :);
rBand2 = T2(T2.series == "Real kart +/-15% band" & T2.x_m <= 114, :);

figure(2);
hold on; grid on;

% Plot lines (change colors/styles here)
plot(dyn2.x_m,  dyn2.speed_mps,  'b',   'LineWidth', 1.5);
plot(eb2.x_m,   eb2.speed_mps,   'g--', 'LineWidth', 1.5);
plot(real2.x_m, real2.speed_mps, 'k',   'LineWidth', 1.5);

% +/- 15% bounds
% plot(dBand2.x_m, dBand2.lower_speed_mps, 'b:');
% plot(dBand2.x_m, dBand2.upper_speed_mps, 'b:', 'HandleVisibility', 'off');
% 
% plot(rBand2.x_m, rBand2.lower_speed_mps, 'g:');
% plot(rBand2.x_m, rBand2.upper_speed_mps, 'g:', 'HandleVisibility', 'off');

xlim([0 113]);
xlabel('Distance (m)');
ylabel('Speed (m/s)');
title('Course Speed Comparison');
legend('Dynamic', 'Energy balance', 'Real kart speed', 'Location', 'best');
hold off;