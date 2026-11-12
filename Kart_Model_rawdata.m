clear; clc; close all;

%% 1. Raw GPS Data
mean_profile_data = [
    0  1.13686837722e-13;    1  -0.141326132126;    2  -0.103689499265;
    3  0.0209335309806;     4  0.110879711346;     5  0.113225563153;
    6  0.323557598781;      7  0.454000107093;     8  0.614167411632;
    9  0.791423546685;     10  1.03530777571;     11  1.3337362667;
   12  1.58295130124;      13  1.77886351526;     14  1.95303156407;
   15  2.01729080074;      16  2.16702573292;     17  2.32483076454;
   18  2.54211248135;      19  2.70807833094;     20  2.82502011356;
   21  2.90586649131;      22  2.94969012622;     23  3.00993575444;
   24  3.08760530618;      25  3.1642922574;      26  3.23017850077;
   27  3.27384955897;      28  3.31033030486;     29  3.34134672753;
   30  3.36872715794;      31  3.39137341993;     32  3.48289699338;
   33  3.65921275247;      34  3.87399623129;     35  3.99175376111;
   36  4.05623520409;      37  4.07147749062;     38  4.07161006656;
   39  4.06583122759;      40  4.11487712939;     41  4.23070391746;
   42  4.36857465945;      43  4.48395054869;     44  4.51317575575;
   45  4.56252118869;      46  4.61777033729;     47  4.62830568873;
   48  4.67299211188;      49  4.73435515568;     50  4.78245809414;
   51  4.88738795105;      52  5.01805621771;     53  5.16053878987;
   54  5.30091354946;      55  5.3740642101;      56  5.41445046666;
   57  5.44349310442;      58  5.42044683105;     59  5.39161121072;
   60  5.3825079793;       61  5.3719676412;      62  5.32373635107;
   63  5.29532186711;      64  5.30744263011;     65  5.30870585555;
   66  5.43483038578;      67  5.52529153366;     68  5.60005147792;
   69  5.64531441882;      70  5.69731868662;     71  5.74992647024;
   72  5.79265205967;      73  5.94824485463;     74  6.14065308553;
   75  6.25799148261;      76  6.4416385292;      77  6.64203701104;
   78  6.83292086709;      79  6.96634270603;     80  7.0186585725;
   81  7.00789447575;      82  7.04312237346;     83  7.03460290447;
   84  7.07165321739;      85  7.17163670328;     86  7.27277466224;
   87  7.32985483083;      88  7.32320686882;     89  7.27697118214;
   90  7.2326332833;       91  7.18696960709;     92  7.13648355517;
   93  7.10886243252;      94  7.13793188229;     95  7.17000260299;
   96  7.17535699728;      97  7.12381863051;     98  6.9948912124;
   99  6.78786036072;     100  6.60331476616;    101  6.49155197866;
  102  6.52942608396;     103  6.67259118907;    104  6.88550502092;
  105  7.17598370373;     106  7.30317913775;    107  7.37155962213;
  108  7.41934251384;     109  7.45684616246;    110  7.52054760107;
  111  7.66377576071;     112  7.90762862967;    113  8.13065223529;
  114  8.33215324193;     115  8.4917131227;     116  8.75617491295;
  117  9.04350389375;     118  9.33948005532;    119  9.56831524898;
  119.979435012  9.82015802558;
];

% Configure course: start at top (x=0, y=9.82 m) and finish at bottom (x=119.98 m, y=0 m)
x_data = mean_profile_data(:, 1);
y_data = flipud(mean_profile_data(:, 2)); 

%% 2. Fit Cubic Smoothing Spline
% p = 0.999 provides smooth 1st and 2nd derivatives without numerical ripple
track_spline = csaps(x_data, y_data, 0.999);
d_spline     = fnder(track_spline, 1);
dd_spline    = fnder(track_spline, 2);

f   = @(x) fnval(track_spline, x);
df  = @(x) fnval(d_spline, x);
ddf = @(x) fnval(dd_spline, x);

%% 3. Vehicle Parameters
M     = 136.078;                   % Total mass (car + driver) [kg]
R     = 0.127;                 % Wheel radius [m] (~10 in diameter)
Iw    = 0.015;                 % Polar inertia per wheel [kg*m^2]
Meff  = M + 4 * (Iw / R^2);     % Effective mass [kg]
g     = 9.81;                  % Gravity [m/s^2]
Crr   = 0.03;                 % Rolling resistance coefficient [-]

%% 4. ODE Function
v_smooth = 0.05; % Smoothing parameter to prevent zero-velocity singularity

odefun = @(t, z) [
    z(2);
    (-M * g * df(z(1)) ...
     - Meff * df(z(1)) * ddf(z(1)) * (z(2)^2) ...
     - Crr * M * g * tanh(z(2) / v_smooth) / sqrt(1 + (df(z(1))^2))) / ...
    (Meff * (1 + (df(z(1))^2)))
];

%% 5. Solve ODE
x_start = x_data(1);
x_end   = x_data(end);
z0      = [x_start; 0];
tspan   = [0, 40];

% Stop solver exactly when passing the finish line
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, ...
                 'Events', @(t, z) finishLineEvent(t, z, x_end));

[t, z] = ode45(odefun, tspan, z0, options);

%% 6. Post-Process Results
x     = z(:, 1);
x_dot = z(:, 2);
y     = f(x);
y_dot = df(x) .* x_dot;
v     = sqrt(x_dot.^2 + y_dot.^2); % Total tangential speed [m/s]
v_mph = v * 2.23694;

fprintf('========================================\n');
fprintf('Finish Line Time: %.2f seconds\n', t(end));
fprintf('Finish Line Speed: %.2f m/s (%.1f mph)\n', v(end), v_mph(end));
fprintf('Max Speed Reached: %.2f m/s (%.1f mph)\n', max(v), max(v_mph));
fprintf('========================================\n');

%% 7. Diagnostic Plots
figure('Name', 'Soapbox Run Performance', 'Color', 'w', 'Position', [100, 100, 1000, 700]);

subplot(2, 2, 1);
plot(x_data, y_data, 'k.', 'MarkerSize', 8, 'DisplayName', 'GPS Data'); hold on;
x_fine = linspace(x_start, x_end, 500);
plot(x_fine, f(x_fine), 'b-', 'LineWidth', 2, 'DisplayName', 'Spline Fit f(x)');
grid on; xlabel('Distance x [m]'); ylabel('Elevation y [m]');
title('Track Elevation Profile'); legend('Location', 'best');

subplot(2, 2, 2);
plot(x_fine, df(x_fine), 'r-', 'LineWidth', 1.5);
yline(0, 'k--');
grid on; xlabel('Distance x [m]'); ylabel("Slope f'(x)");
title("Track Slope Profile");

subplot(2, 2, 3);
plot(t, v_mph, 'Color', [0.13, 0.55, 0.13], 'LineWidth', 2);
grid on; xlabel('Time [s]'); ylabel('Speed [mph]');
title('Velocity vs. Time');

subplot(2, 2, 4);
plot(x, v_mph, 'm-', 'LineWidth', 2);
grid on; xlabel('Distance x [m]'); ylabel('Speed [mph]');
title('Velocity vs. Track Distance');

%% Event Function
function [value, isterminal, direction] = finishLineEvent(~, z, x_finish)
    value      = z(1) - x_finish;
    isterminal = 1;
    direction  = 1;
end