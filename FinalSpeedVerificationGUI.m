function KartVelocitySimulation_v51
% KARTVELOCITYSIMULATION_V51 Extended hill with 50% slower angle variation.
% Save as KartVelocitySimulation_v51.m and run KartVelocitySimulation_v51.
% Requires MATLAB with uifigure/uigridlayout support; no extra toolboxes.
%
% Changelog
% v51: Use linear b*v resistance in both regions and both solvers.
% v51: Based on the uploaded v36 file. End at section x=25 m.
%      thetaNew(x)=thetaStart+0.5*(thetaOriginal(x)-thetaStart).
%      d(thetaNew)/dx is exactly half d(thetaOriginal)/dx.
%      Integrate the new slope for consistent height in both solvers.
% v36: Based on the supplied v34(1) file. Use quadratic resistance b*v^2
%      in both friction regions and both solvers. b has units N*s^2/m^2.
%      Dissipated power is b*v^3; update wheel coefficients and UI units.
%      Preserve supplied parameters, hill, offsets, and acceleration equation.
% v34: Shift simulated comparison curves and the dynamic band left 0.4 m.
%      Trim comparison samples before x=0, including selected CSV exports.
%      Preserve physical trajectories, hill profile, and GPS coordinates.
% v33: Replace the speed fit and both raw recordings with the averaged CSV
%      GPS Speed samples (0..30 m), embedded for standalone use. Keep the
%      assumed +/-15% GPS band, visibility controls, and selected exports.
%      Compare models with GPS at matched comparison coordinates.
% v32: Shift Raw 1 left by 5 m and Raw 2 left by 2.52 m, then discard
%      samples at negative shifted coordinates. Draw raw curves without markers.
% v31: Restore the original vertical position of the real-speed polynomial.
%      Restore selectable raw GPS curves and +/-15% bands with selected export.
% v30: Add b1 before the hill-coordinate switch at x=5 m and b2 after it.
%      Shift only speed-comparison curves and bands right by 0.8 m.
%      Hill coordinates and simulation integration coordinates are unchanged.
% v29: Restore separate mu1 and mu2 inputs before/after x=5 m.
%      Display both simulated speeds with a -0.15 m/s vertical offset.
%      Shift the dynamic band with its curve; retain physical solver speeds.
% v28: Move the friction switch to x=5 m. Export checked comparison
%      curves and bands as CSV data or the selected plot as PNG/PDF.
% v27: Remove the raw GPS recordings and add a band around the real-speed fit.
%      Use one Coulomb coefficient in both regions. Shift the entire hill
%      left by 0.4 m, including the endpoints and friction switch.
% v26: Restore the endpoint to 20.93 m and set the speed fit to zero at xStart.
%      Add separately selectable +/-15% bands for both raw GPS speed curves.
% v25: Extend the simulation endpoint to displayed x=40 m.
%      Continue the hill polynomial beyond the previous x=20.93 m endpoint.
%      Shift the speed fit vertically to 0.2 m/s at xStart=-0.4 m.
% v24: Add comparison visibility checkboxes and an assumed +/-15% speed band.
%      Extend the hill left by 0.4 m to displayed x=-0.4 m.
% v23: Add adjustable Coulomb coefficient mu2 after the friction switch.
%      Second-region resistance is b*v + mu2*m*g*cos(theta).
%      Include both Coulomb regions in energy losses and CSV exports.
% v22: Begin at source x=14.07 m, previous displayed x=2.65 m.
%      This is just beyond the mu=0.03 acceleration threshold (14.06313 m).
%      Rebase to displayed 0..20.93 m, with the friction switch at 7.35 m.
% v21: Offset simulated coordinates by -10 m to align the start with GPS x=0.
%      Use displayed 0..23.58 m for the same hill source interval 11.42..35 m.
%      The friction transition is now at x=10 m. Initial speed is adjustable.
% v20: Restrict both solvers and plots to displayed x=10..33.58 m.
%      Initial velocity and zero relative elevation apply at x=10 m.
% v19: Rebase the previous displayed -10 m to zero. Start both solvers
%      at zero using the original hill interval 1.42..35 m (33.58 m long).
%      Add adjustable initial speed. Move the previous 10 m friction
%      boundary to 20 m on the rebased axis.
% v18: Start at displayed x=-4 m (original x=6 m) with the fitted speed there.
% v17: Start both solvers from rest at displayed x=-10 m (original x=0 m).
%      Rebase relative elevation and elapsed distance at the new start.
% Piecewise friction: use Coulomb resistance before displayed x=10 m
% and viscous resistance from x=10 m onward in both solvers.
% v11: Embed both supplied GPS recordings and overlay unfiltered recorded
%      speeds as connected lines. Add a full-size Speed comparison tab.
%      GPS distance is cumulative horizontal great-circle distance from
%      each recording's first fix; displayed x=GPS distance-10 m.
%      This origin assumes the first fix is the original course origin.
%      CSV raw-speed columns are linearly interpolated at matching x.
% v10: Start at displayed 4 m (original 14 m), finish at displayed 25 m.
%      Fix initial speed to the reference value there: 1.855699280 m/s.
% v9: Remove initial-speed control; fix v0=1.974809375 m/s in both models.
% v8: Restrict BOTH solvers and reference curves to current x=5..25 m
%     (original x=15..35 m). Keep the displayed coordinate at 5..25 m.
%     Reset time/path distance/losses at x=5; match default initial speed
%     to the real-speed fit at original 15 m. Height is zero at this start.
% v7: Add adjustable total wheel viscous resistance Fviscous=cViscous*v.
%     Apply the force in dynamics and its work integral in energy balance.
%     Export viscous force/loss and plot cumulative viscous energy loss.
%     Default cViscous=0 preserves v6; positive values enable damping.
% v6: Simplify controls to mass, rolling resistance, initial speed, and time
%     limit. Fix course/wheel/air settings and neglect aerodynamic drag in
%     both models (CdA=0). Keep zero drag columns in CSVs for compatibility.
% v5: Restrict both curves to original x=10..35 m; rebase horizontal x to
%     0..25 m and elevation to zero at the new start. Preserve actual speed
%     values; default initial speed equals the reference fit at original 10 m.
%     Export original and local coordinates; correct energy-label newlines.
% v4: Independent work-energy integration in horizontal x, three speed
%     overlays, energy accounting, model agreement, and energy CSV export.
% v3: Overlay supplied real-speed polynomial against horizontal x; export
%     reference speed and simulated-minus-reference error at matched x.
%     The supplied fit is shown unchanged, including its negative intercept.
% v2: Polynomial elevation profile, local slope, horizontal x bounds,
%     along-path distance integration, hill plots, and profile CSV columns.
% v1: Original constant-slope GUI (preserved as KartVelocitySimulation.m).
%
% Assumptions: x and elevation are in meters; travel toward increasing x;
% rolling without slip; negligible aerodynamic drag.
% Before the switch, F = b1*v + mu1*m*g*cos(theta).
% From the switch onward, F = b2*v + mu2*m*g*cos(theta).
% mu1 and mu2 are stored as Crr and Crr2. b1 and b2 are stored as
% cViscous1 and cViscous (N*s/m), representing all wheels together.
% Linear resistance power is b*v^2. For each wheel, torque = k*omega,
% where k = b*radius^2/nWheels (N*m*s). Coefficients describe all wheels.
% Numeric b defaults are preserved from the supplied file; units have changed.
% Curvature-induced normal-load changes, airborne motion, and rollback are
% not modeled. Source x=13.67..39.47 m maps to section x=-0.8..25 m.
% The original polynomial is extrapolated beyond its previous source end 35 m.
% Mass and friction defaults are examples, not fitted measurements.
% Hlocal(x) is the integral of the modified slope from the original start.
% Initial slope and elevation datum are preserved; the total height drop changes.
% GPS Speed uses the CSV x_m and speed_mps values without added offsets.
% GPS shading is the assumed +/-15% speed band, not the CSV run envelope.
% Initial speed is set by the user and applies to both solvers.

%% Defaults and input definitions
% key, label, default, limits, tooltip
spec = {
    'mass',       'Kart + rider mass (kg)',        117,    [0.001 Inf], 'Total mass INCLUDING all wheels.';
    'Crr',        'Coulomb friction 1, μ₁',       0.02,   [0 Inf],     'Dimensionless coefficient used before hill x = 11.6 m. FCoulomb = mu1*m*g*cos(theta).';
    'Crr2',       'Coulomb friction 2, μ₂',       0.03,   [0 Inf],     'Dimensionless coefficient used from hill x = 11.6 m onward. F = b2*v + mu2*m*g*cos(theta).';
    'cViscous1',  'Linear b₁ (N*s/m)',          15,      [0 Inf],     'Linear resistance before hill x = 11.6 m. Fviscous = b1*v.';
    'cViscous',   'Linear b₂ (N*s/m)',          15,      [0 Inf],     'Linear resistance from hill x = 11.6 m onward. Fviscous = b2*v.';
    'v0',         'Initial velocity (m/s)',       0,      [0 Inf],     'Speed along the hill at displayed x = -0.8 m, used by both solvers.';
    'tMax',       'Dynamic time limit (s)',             120,    [0.01 3600], 'Time limit for the dynamic solver. The energy method independently integrates the selected x interval until it finishes or stops.'
    };
% Fixed model settings (not GUI controls).
fixed = struct('xStart',-0.8,'xEnd',25,'angleVariationScale',0.5,'xFrictionSwitch',11.6,'sourceXOffset',14.47, ...
    'CdA',0,'rho',1.225,'nWheels',4,'radius',0.127,'Iwheel',0.0175, ...
    'speedErrorFraction',0.15,'simulatedSpeedOffset_mps',-0.15, ...
    'comparisonXOffset_m',0.4); % v33 offset 0.8 minus a 0.4 m left shift.
nInputs = size(spec,1);
fields = struct();
lastResults = table();
lastParameters = struct();
lastEnergy = table();
gpsSpeed = embeddedGPSSpeed();
busy = false;
blue = [0.08 0.38 0.68];

%% Window and parameter panel
fig = uifigure('Name','Kart Velocity Simulation v51 (-0.8 to 25 m)', ...
    'Position',[80 80 1250 780], 'Color',[0.96 0.97 0.98]);
outer = uigridlayout(fig,[1 2]);
outer.ColumnWidth = {335,'1x'};
outer.Padding = [14 14 14 14];
outer.ColumnSpacing = 14;

panel = uipanel(outer,'Title','Parameters','FontWeight','bold');
controls = uigridlayout(panel,[nInputs+7 2]);
controls.ColumnWidth = {'1x',100};
controls.RowHeight = [repmat({36},1,nInputs), {28,38,34,34,130,34,'1x'}];
controls.Padding = [12 12 12 12];
controls.RowSpacing = 7;
controls.Scrollable = 'on';
for i = 1:size(spec,1)
    label = uilabel(controls,'Text',spec{i,2},'FontSize',12);
    label.Layout.Row = i; label.Layout.Column = 1;
    f = uieditfield(controls,'numeric','Value',spec{i,3}, ...
        'Limits',spec{i,4},'Tooltip',spec{i,5}, ...
        'ValueChangedFcn',@parametersChanged);
    f.Layout.Row = i; f.Layout.Column = 2;
    fields.(spec{i,1}) = f;
end

autoRun = uicheckbox(controls,'Text','Update plots after each edit', ...
    'Value',true,'ValueChangedFcn',@parametersChanged);
autoRun.Layout.Row = nInputs+1; autoRun.Layout.Column = [1 2];
runButton = uibutton(controls,'Text','Run simulation','FontWeight','bold', ...
    'BackgroundColor',blue,'FontColor',[1 1 1], 'ButtonPushedFcn',@runSimulation);
runButton.Layout.Row = nInputs+2; runButton.Layout.Column = [1 2];
resetButton = uibutton(controls,'Text','Reset example values','ButtonPushedFcn',@resetInputs);
resetButton.Layout.Row = nInputs+3; resetButton.Layout.Column = [1 2];
exportButton = uibutton(controls,'Text','Export results to CSV', ...
    'Enable','off','ButtonPushedFcn',@exportResults);
exportButton.Layout.Row = nInputs+4; exportButton.Layout.Column = [1 2];
note = uilabel(controls,'Text',sprintf(['Course x = -0.8 to 25 m; angle variation reduced 50%.\n' ...
    'Before hill x = 11.6 m: F = b1*v + mu1*m*g*cos(theta).\n' ...
    'From hill x = 11.6 m: F = b2*v + mu2*m*g*cos(theta).\n' ...
    'Initial velocity is set above.\n' ...
    'Aerodynamic drag is neglected.']), ...
    'WordWrap','on','FontColor',[0.35 0.39 0.44],'VerticalAlignment','top');
note.Layout.Row = nInputs+5; note.Layout.Column = [1 2];
energyExportButton = uibutton(controls,'Text','Export full energy CSV', ...
    'Enable','off','ButtonPushedFcn',@exportEnergy);
energyExportButton.Layout.Row = nInputs+6; energyExportButton.Layout.Column = [1 2];

%% Results area
right = uigridlayout(outer,[3 1]);
right.RowHeight = {82,'1x',55};
right.Padding = [0 0 0 0];
metrics = uigridlayout(right,[1 3]);
metrics.Padding = [0 0 0 0];
metricSpeed = uilabel(metrics,'Text','Final speed','FontSize',16,'FontWeight','bold');
metricTime = uilabel(metrics,'Text','Elapsed time','FontSize',16,'FontWeight','bold');
metricMass = uilabel(metrics,'Text','Effective mass','FontSize',16,'FontWeight','bold');
tabs = uitabgroup(right);
comparisonTab = uitab(tabs,'Title','Speed comparison');
comparisonGrid = uigridlayout(comparisonTab,[3 1]);
comparisonGrid.RowHeight = {100,'1x',76};
visibilityPanel = uipanel(comparisonGrid,'Title','Show curves');
visibilityPanel.Layout.Row = 1;
visibilityGrid = uigridlayout(visibilityPanel,[2 3]);
visibilityGrid.RowHeight = {24,24};
visibilityGrid.RowSpacing = 5;
visibilityGrid.Padding = [8 8 8 8];
visibilityGrid.ColumnWidth = {'1x','1x','1x'};
comparisonLabels = {'Dynamic','Energy balance','GPS Speed', ...
    'Dynamic +/-15% band','GPS +/-15% band'};
comparisonChecks = gobjects(1,numel(comparisonLabels));
for k = 1:numel(comparisonLabels)
    comparisonChecks(k) = uicheckbox(visibilityGrid,'Text',comparisonLabels{k}, ...
        'Value',true,'ValueChangedFcn',@updateComparisonVisibility);
    comparisonChecks(k).Layout.Row = ceil(k/3);
    comparisonChecks(k).Layout.Column = mod(k-1,3)+1;
end
comparisonFooter = uigridlayout(comparisonGrid,[1 2]);
comparisonFooter.Layout.Row = 3;
comparisonFooter.ColumnWidth = {'1x',160};
comparisonFooter.Padding = [0 0 0 0];
comparisonExportButton = uibutton(comparisonFooter,'Text','Export selected...', ...
    'Enable','off','ButtonPushedFcn',@exportComparison, ...
    'Tooltip','Export all samples of checked curves/bands as CSV, or the current plot view as PNG/PDF.');
comparisonExportButton.Layout.Row = 1;
comparisonExportButton.Layout.Column = 2;
axComparison = uiaxes(comparisonGrid);
axComparison.Layout.Row = 2;
bandNote = uilabel(comparisonFooter,'Text', ...
    sprintf(['Shading shows the assumed +/-15%% speed bands.\n' ...
    'Simulation: -0.15 m/s; shifted left 0.4 m; cropped before x = 0.\n' ...
    'GPS uses the CSV coordinates directly (0-30 m), without added offsets.']), ...
    'WordWrap','on','FontSize',11);
bandNote.Layout.Row = 1;
bandNote.Layout.Column = 1;
motionTab = uitab(tabs,'Title','Motion');
hillTab = uitab(tabs,'Title','Hill profile');
energyTab = uitab(tabs,'Title','Energy balance');
energyGrid = uigridlayout(energyTab,[4 1]);
energyGrid.RowHeight = {66,'1x','1x',60};
energyNote = uilabel(energyGrid,'Text',sprintf([ ...
    'K = K0 + m*g*(h0-h) - WCoulomb - Wviscous;  K = 0.5*mEffective*v^2\n' ...
    'Before hill x = 11.6 m: mu1 and b1. From hill x = 11.6 m: mu2 and b2.\n' ...
    'Course: -0.8 to 25 m; downhill angle rate reduced 50%. Time limit applies only to dynamics.']), ...
    'WordWrap','on','FontSize',12);
axEnergy = uiaxes(energyGrid);
axAgreement = uiaxes(energyGrid);
energySummary = uilabel(energyGrid,'Text','Ready.','WordWrap','on');
plotGrid = uigridlayout(motionTab,[2 2]);
plotGrid.Padding = [0 0 0 0];
axSpeed = uiaxes(plotGrid); axSpeed.Layout.Row = 1; axSpeed.Layout.Column = 1;
axDistanceSpeed = uiaxes(plotGrid); axDistanceSpeed.Layout.Row = 1; axDistanceSpeed.Layout.Column = 2;
axAccel = uiaxes(plotGrid); axAccel.Layout.Row = 2; axAccel.Layout.Column = 1;
axDistance = uiaxes(plotGrid); axDistance.Layout.Row = 2; axDistance.Layout.Column = 2;
hillGrid = uigridlayout(hillTab,[3 1]);
hillGrid.RowHeight = {48,'1x','1x'};
formulaLabel = uilabel(hillGrid,'Text', ...
    'theta(x) = thetaStart + 0.5*(thetaOriginal(x)-thetaStart); height integrates the new slope. Course: -0.8 to 25 m.', ...
    'WordWrap','on','FontSize',12);
axHill = uiaxes(hillGrid);
axGrade = uiaxes(hillGrid);
statusLabel = uilabel(right,'Text','Ready.','WordWrap','on', ...
    'FontColor',[0.30 0.34 0.40]);
runSimulation();

%% Callbacks
    function updateComparisonVisibility(varargin)
        handles = axComparison.UserData;
        if isempty(handles) || ~all(isgraphics(handles)), return; end
        selected = false(size(handles));
        for j = 1:numel(handles)
            selected(j) = comparisonChecks(j).Value;
            if selected(j)
                handles(j).Visible = 'on';
            else
                handles(j).Visible = 'off';
            end
        end
        if any(selected)
            legend(axComparison,handles(selected),comparisonLabels(selected), ...
                'Location','southoutside','NumColumns',3,'FontSize',11, ...
                'AutoUpdate','off');
        else
            legend(axComparison,'off');
        end
        comparisonExportButton.Enable = 'off';
        if any(selected) && strcmp(exportButton.Enable,'on')
            comparisonExportButton.Enable = 'on';
        end
    end

    function parametersChanged(~,~)
        exportButton.Enable = 'off';
        energyExportButton.Enable = 'off';
        comparisonExportButton.Enable = 'off';
        if autoRun.Value
            runSimulation();
        else
            statusLabel.Text = 'Inputs changed. Press Run simulation to update these plots.';
        end
    end

    function resetInputs(~,~)
        for j = 1:size(spec,1)
            fields.(spec{j,1}).Value = spec{j,3};
        end
        runSimulation();
    end

    function runSimulation(varargin)
        if busy, return; end
        busy = true;
        runButton.Enable = 'off';
        exportButton.Enable = 'off';
        energyExportButton.Enable = 'off';
        comparisonExportButton.Enable = 'off';
        cleanup = onCleanup(@finishRun); %#ok<NASGU>
        statusLabel.Text = 'Calculating...';
        drawnow limitrate nocallbacks
        try
            p = fixed;
            for j = 1:size(spec,1)
                p.(spec{j,1}) = fields.(spec{j,1}).Value;
            end
            vals = struct2arrayLocal(p);
            if any(~isfinite(vals))
                error('All inputs must be finite numbers.');
            end
            [result,p,message] = simulateKart(p);
            [energy,energyMessage] = simulateEnergy(p);
            % Compare at the same x, not at the same row or solver time.
            result.energy_kinetic_J = sampleAtX(energy.x_m,energy.kinetic_J,result.x_m);
            result.energy_speed_mps = sqrt(2*max(result.energy_kinetic_J,0)/p.mEffective);
            result.energy_speed_mps(~isfinite(result.energy_kinetic_J)) = NaN;
            result.dynamic_minus_energy_mps = result.speed_mps-result.energy_speed_mps;
            % Apply plotting offsets after solving the physical motion.
            result.comparison_x_m = result.x_m+p.comparisonXOffset_m;
            energy.comparison_x_m = energy.x_m+p.comparisonXOffset_m;
            % Linear interpolation matches the straight segments of the GPS curve.
            % NaN outside the supplied CSV range; never extrapolate GPS speeds.
            result.gps_speed_mps = interp1(gpsSpeed.x_m,gpsSpeed.speed_mps, ...
                result.comparison_x_m,'linear',NaN);
            energy.gps_speed_mps = interp1(gpsSpeed.x_m,gpsSpeed.speed_mps, ...
                energy.comparison_x_m,'linear',NaN);
            result.speed_error_mps = result.speed_mps-result.gps_speed_mps;
            result.energy_minus_gps_mps = result.energy_speed_mps-result.gps_speed_mps;
            energy.energy_minus_gps_mps = energy.speed_mps-energy.gps_speed_mps;
            result.displayed_speed_mps = result.speed_mps+p.simulatedSpeedOffset_mps;
            energy.displayed_speed_mps = energy.speed_mps+p.simulatedSpeedOffset_mps;
            result.displayed_speed_error_mps = result.displayed_speed_mps-result.gps_speed_mps;
            energy.displayed_speed_error_mps = energy.displayed_speed_mps-energy.gps_speed_mps;
            result.dynamic_speed_lower_mps = (1-p.speedErrorFraction)*result.speed_mps+p.simulatedSpeedOffset_mps;
            result.dynamic_speed_upper_mps = (1+p.speedErrorFraction)*result.speed_mps+p.simulatedSpeedOffset_mps;
            result.gps_speed_lower_mps = result.gps_speed_mps-p.speedErrorFraction*abs(result.gps_speed_mps);
            result.gps_speed_upper_mps = result.gps_speed_mps+p.speedErrorFraction*abs(result.gps_speed_mps);
            energy.gps_speed_lower_mps = energy.gps_speed_mps-p.speedErrorFraction*abs(energy.gps_speed_mps);
            energy.gps_speed_upper_mps = energy.gps_speed_mps+p.speedErrorFraction*abs(energy.gps_speed_mps);
            lastEnergy = energy;
            lastResults = result;
            lastParameters = p;
            fig.UserData = struct('Results',result,'Energy',energy,'Parameters',p,'GPSSpeed',gpsSpeed);
            plotResult(axSpeed,result.time_s,result.displayed_speed_mps, ...
                'Dynamic velocity over time','Time (s)','Speed (m/s)');
            drawComparison(axDistanceSpeed,result,energy,p,false);
            drawComparison(axComparison,result,energy,p,true);
            plotResult(axAccel,result.time_s,result.acceleration_mps2, ...
                'Acceleration over time','Time (s)','Acceleration (m/s^2)');
            plotResult(axDistance,result.time_s,result.distance_m, ...
                'Distance over time','Time (s)','Distance along hill (m)');
            plotResult(axEnergy,energy.x_m,energy.kinetic_J, ...
                'Work-energy accounting','Section coordinate x (m)','Energy (J)');
            hold(axEnergy,'on');
            plot(axEnergy,energy.x_m,energy.initial_kinetic_J+energy.gravity_work_J, ...
                ':','Color',[0.45 0.30 0.70],'LineWidth',1.6);
            plot(axEnergy,energy.x_m,energy.rolling_loss_J+energy.viscous_loss_J, ...
                '--','Color',[0.90 0.35 0.08],'LineWidth',1.6);
            plot(axEnergy,energy.x_m,energy.viscous_loss_J,'-.', ...
                'Color',[0.65 0.18 0.48],'LineWidth',1.6);
            hold(axEnergy,'off');
            legend(axEnergy,{'Kinetic energy','Initial kinetic + gravity work', ...
                'Coulomb + viscous losses','Viscous loss only'},'Location','best');
            plotResult(axAgreement,result.x_m,result.dynamic_minus_energy_mps, ...
                'Dynamic speed minus energy speed','Section coordinate x (m)','Speed difference (m/s)');
            difference = result.dynamic_minus_energy_mps;
            difference = difference(isfinite(difference));
            if isempty(difference), maxDifference = NaN; else, maxDifference = max(abs(difference)); end
            energySummary.Text = sprintf(['%s\n' ...
                'Max speed difference on shared x: %.3g m/s | Max energy residual: %.3g J'], ...
                energyMessage,maxDifference,max(abs(energy.balance_residual_J)));
            xHill = linspace(p.xStart,p.xEnd,600);
            elevationHill = hillHeight(xHill,p);
            plotResult(axHill,xHill,elevationHill, ...
                'Hill elevation relative to new start','Section coordinate x (m)','Relative elevation (m)');
            hold(axHill,'on');
            plot(axHill,result.x_m,result.elevation_m,'Color',[0.95 0.45 0.12],'LineWidth',2);
            plot(axHill,result.x_m(end),result.elevation_m(end),'o', ...
                'MarkerFaceColor',[0.95 0.45 0.12],'MarkerEdgeColor','w','MarkerSize',8);
            hold(axHill,'off');
            legend(axHill,{'Hill profile','Simulated path','Final position'},'Location','best');
            plotResult(axGrade,xHill,-atan(hillSlope(xHill,p))*180/pi, ...
                'Local slope (positive downhill)','Section coordinate x (m)','Downhill angle (deg)');
            xlim(axEnergy,[p.xStart p.xEnd]);
            xlim(axAgreement,[p.xStart p.xEnd]);
            xlim(axHill,[p.xStart p.xEnd]);
            xlim(axGrade,[p.xStart p.xEnd]);
            metricSpeed.Text = sprintf('Displayed final speed\n%.2f m/s | %.2f mph', ...
                result.displayed_speed_mps(end),result.displayed_speed_mps(end)*2.236936292);
            metricTime.Text = sprintf('Elapsed time\n%.2f s | %.2f m', ...
                result.time_s(end),result.distance_m(end));
            metricMass.Text = sprintf('Effective mass\n%.2f kg (+%.2f)', ...
                p.mEffective,p.mEffective-p.mass);
            statusLabel.Text = message;
            exportButton.Enable = 'on';
            energyExportButton.Enable = 'on';
            updateComparisonVisibility();
        catch err
            statusLabel.Text = 'Simulation failed. Check the inputs; any existing plots show the previous run.';
            uialert(fig,err.message,'Simulation error');
        end
    end

    function finishRun()
        busy = false;
        if isvalid(fig), runButton.Enable = 'on'; end
    end

    function plotResult(ax,x,y,heading,xtext,ytext)
        plot(ax,x,y,'Color',blue,'LineWidth',2);
        grid(ax,'on');
        title(ax,heading); xlabel(ax,xtext); ylabel(ax,ytext);
        ax.FontSize = 11;
    end

    function drawComparison(ax,result,energy,p,large)
        % Crop plotting copies only. Full solver tables retain physical history.
        [dynamicX,dynamicValues] = trimComparisonSamples(result.comparison_x_m, ...
            [result.displayed_speed_mps result.dynamic_speed_lower_mps result.dynamic_speed_upper_mps]);
        [energyX,energySpeed] = trimComparisonSamples(energy.comparison_x_m,energy.displayed_speed_mps);
        cla(ax);
        hold(ax,'on');
        gpsColor = [0.90 0.35 0.08];
        gpsX = gpsSpeed.x_m;
        gpsV = gpsSpeed.speed_mps;
        gpsLower = gpsV-p.speedErrorFraction*abs(gpsV);
        gpsUpper = gpsV+p.speedErrorFraction*abs(gpsV);
        if large
            hBand = fill(ax,[dynamicX;flipud(dynamicX)], ...
                [dynamicValues(:,2);flipud(dynamicValues(:,3))], ...
                blue,'FaceAlpha',0.18,'EdgeColor','none', ...
                'HitTest','off','PickableParts','none');
            hGPSBand = fill(ax,[gpsX;flipud(gpsX)], ...
                [gpsLower;flipud(gpsUpper)],gpsColor, ...
                'FaceAlpha',0.18,'EdgeColor','none', ...
                'HitTest','off','PickableParts','none');
        end
        hDynamic = plot(ax,dynamicX,dynamicValues(:,1),'-', ...
            'Color',blue,'LineWidth',2.2,'Marker','none');
        hEnergy = plot(ax,energyX,energySpeed,'--', ...
            'Color',[0.05 0.60 0.36],'LineWidth',2.0,'Marker','none');
        hGPS = plot(ax,gpsX,gpsV,'-', ...
            'Color',gpsColor,'LineWidth',2.2,'Marker','none');
        hold(ax,'off');
        grid(ax,'on');
        % Show all GPS samples without extending the simulated hill interval.
        xlim(ax,[min(0,gpsX(1)), ...
            max(p.xEnd+p.comparisonXOffset_m,gpsX(end))]);
        visibleValues = [dynamicValues(:,1);energySpeed;gpsV];
        if large
            visibleValues = [visibleValues;dynamicValues(:,2); ...
                dynamicValues(:,3);gpsLower;gpsUpper];
        end
        visibleValues = visibleValues(isfinite(visibleValues));
        top = max(visibleValues); bottom = min([0;visibleValues]);
        ylim(ax,[bottom-0.02*max(1,top-bottom),top+0.08*max(1,top-bottom)]);
        title(ax,'Simulated and GPS speed');
        xlabel(ax,'Comparison coordinate (m)'); ylabel(ax,'Speed (m/s)');
        labels = {'Dynamic','Energy balance','GPS Speed'};
        if large
            ax.FontSize = 13;
            ax.UserData = [hDynamic hEnergy hGPS hBand hGPSBand];
            updateComparisonVisibility();
        else
            ax.FontSize = 10;
            legend(ax,[hDynamic hEnergy hGPS],labels,'Location','best','FontSize',8);
        end
    end

    function exportComparison(~,~)
        if ~strcmp(comparisonExportButton.Enable,'on'), return; end
        filters = {'*.csv','Selected curve and band data (*.csv)'; ...
            '*.png','Selected plot image (*.png)'; ...
            '*.pdf','Selected plot PDF (*.pdf)'};
        [name,folder,formatIndex] = uiputfile(filters,'Export speed comparison', ...
            'KartSpeedComparison_v51.csv');
        if isequal(name,0), return; end
        [~,~,extension] = fileparts(name);
        if isempty(extension)
            extensions = {'.csv','.png','.pdf'};
            extension = extensions{formatIndex};
            name = [name extension];
        end
        destination = fullfile(folder,name);
        try
            switch lower(extension)
                case '.csv'
                    handles = axComparison.UserData;
                    exported = table();
                    for j = 1:numel(handles)
                        if ~comparisonChecks(j).Value, continue; end
                        x = handles(j).XData(:);
                        y = handles(j).YData(:);
                        if isgraphics(handles(j),'line')
                            speed = y;
                            lowerSpeed = nan(size(y));
                            upperSpeed = nan(size(y));
                        else
                            % Each filled polygon stores its lower edge forward
                            % and its upper edge in reverse sample order.
                            count = numel(x)/2;
                            x = x(1:count);
                            speed = nan(size(x));
                            lowerSpeed = y(1:count);
                            upperSpeed = flipud(y(count+1:end));
                        end
                        series = repmat(string(comparisonLabels{j}),numel(x),1);
                        rows = table(series,x,speed,lowerSpeed,upperSpeed, ...
                            'VariableNames',{'series','x_m','speed_mps', ...
                            'lower_speed_mps','upper_speed_mps'});
                        % NaN placeholders represent a simulated curve that
                        % never reached the visible comparison interval.
                        rows = rows(isfinite(rows.x_m),:);
                        exported = [exported;rows]; %#ok<AGROW>
                    end
                    writetable(exported,destination);
                case '.png'
                    exportgraphics(axComparison,destination,'Resolution',300);
                case '.pdf'
                    exportgraphics(axComparison,destination,'ContentType','vector');
                otherwise
                    error('Choose a .csv, .png, or .pdf filename.');
            end
            statusLabel.Text = ['Exported selected comparison to ' destination];
        catch err
            uialert(fig,err.message,'Export error');
        end
    end

    function exportResults(~,~)
        if isempty(lastResults), return; end
        [name,folder] = uiputfile('*.csv','Export simulation results','KartVelocityResults_v51.csv');
        if isequal(name,0), return; end
        try
            % Repeat the run inputs so each CSV is self-contained.
            exported = lastResults;
            keys = fieldnames(lastParameters);
            for j = 1:numel(keys)
                exported.(['input_' keys{j}]) = repmat(lastParameters.(keys{j}),height(exported),1);
            end
            writetable(exported,fullfile(folder,name));
            statusLabel.Text = ['Exported results and run parameters to ' fullfile(folder,name)];
        catch err
            uialert(fig,err.message,'Export error');
        end
    end

    function exportEnergy(~,~)
        if isempty(lastEnergy), return; end
        [name,folder] = uiputfile('*.csv','Export full energy results','KartEnergyResults_v51.csv');
        if isequal(name,0), return; end
        try
            exported = lastEnergy;
            keys = fieldnames(lastParameters);
            for j = 1:numel(keys)
                exported.(['input_' keys{j}]) = repmat(lastParameters.(keys{j}),height(exported),1);
            end
            writetable(exported,fullfile(folder,name));
            statusLabel.Text = ['Exported energy results to ' fullfile(folder,name)];
        catch err
            uialert(fig,err.message,'Export error');
        end
    end
end

%% Dynamics on y = h(x)
function [results,p,status] = simulateKart(p)
% Exact transcription of the supplied polynomial (descending powers).
p.sourceXStart = p.sourceXOffset+p.xStart;
p.sourceXEnd = p.sourceXOffset+p.xEnd;
p.sourceProfileCoefficients = [1.616e-8 -2.653e-6 0.000168 -0.005153 0.07678 -0.5357 1442];
p.sourceElevationAtStart = polyval(p.sourceProfileCoefficients,p.sourceXStart);
p.profileCoefficients = shiftPolynomial(p.sourceProfileCoefficients,p.sourceXOffset);
p.profileCoefficients(end) = p.profileCoefficients(end)-p.sourceElevationAtStart;
p.slopeCoefficients = polyder(p.profileCoefficients); % Original, unmodified slope.
p.startAngle_rad = -atan(polyval(p.slopeCoefficients,p.xStart));
p.mEffective = p.mass + p.nWheels*p.Iwheel/p.radius^2;
p.dragFactor = 0.5*p.rho*p.CdA; % Exactly zero: drag is neglected.
p.linearPerWheel_Nms = p.cViscous*p.radius^2/p.nWheels;
p.linear1PerWheel_Nms = p.cViscous1*p.radius^2/p.nWheels;

% State = [horizontal x; speed along path; distance along path].
% q = sqrt(1 + h'(x)^2), theta = -atan(h'(x)).
% dx/dt = v/q; ds/dt = v.
[initialAcceleration,~,~,~,~] = localForces(p.xStart,p.v0,p,true);
if p.v0 == 0 && initialAcceleration <= 0
    time_s = [0;p.tMax];
    state = [p.xStart 0 0;p.xStart 0 0];
    status = 'Remained at rest: local downhill gravity does not exceed assumed rolling resistance. No rollback is modeled.';
else
    options = odeset('RelTol',1e-8,'AbsTol',[1e-9 1e-10 1e-9],'MaxStep',0.05);
    switchOptions = odeset(options,'Events',@(t,y) stopEvents(t,y,p,p.xFrictionSwitch));
    [time_s,state,~,~,eventIndex] = ode45(@(t,y) kartODE(t,y,p,true), ...
        [0 p.tMax],[p.xStart;p.v0;0],switchOptions);
    reachedEnd = false;
    if any(eventIndex == 1)
        state(end,1) = p.xFrictionSwitch;
        if ~any(eventIndex == 2) && state(end,2) > 0 && time_s(end) < p.tMax
            endOptions = odeset(options,'Events',@(t,y) stopEvents(t,y,p,p.xEnd));
            [viscousTime,viscousState,~,~,eventIndex] = ode45(@(t,y) kartODE(t,y,p,false), ...
                [time_s(end) p.tMax],state(end,:).',endOptions);
            time_s = [time_s;viscousTime(2:end)];
            state = [state;viscousState(2:end,:)];
            reachedEnd = any(eventIndex == 1);
        end
    end
    if reachedEnd
        status = 'Reached the end x coordinate of the selected hill interval.';
    elseif any(eventIndex == 2)
        status = 'Stopped on the profile before reaching end x. No rollback is modeled.';
        state(end,2) = 0;
    else
        status = 'Reached the maximum simulation time before finishing or stopping.';
    end
end

x_m = state(:,1);
distance_m = state(:,3);
elevation_m = hillHeight(x_m,p);
original_x_m = x_m+p.sourceXOffset; % Source coordinate of the hill polynomial.
horizontal_from_start_m = x_m-p.xStart;
original_elevation_m = elevation_m+p.sourceElevationAtStart;
speed_mps = max(state(:,2),0);
speed_mph = speed_mps*2.236936292;
[acceleration_mps2,downhill_angle_deg,gravity_force_N,rolling_force_N,drag_force_N,viscous_force_N] = ...
    localForces(x_m,speed_mps,p,x_m < p.xFrictionSwitch);
acceleration_mps2(speed_mps <= 0 & acceleration_mps2 < 0) = 0;
results = table(time_s,x_m,original_x_m,horizontal_from_start_m,distance_m,elevation_m,original_elevation_m,speed_mps,speed_mph, ...
    acceleration_mps2,downhill_angle_deg,gravity_force_N,rolling_force_N,drag_force_N,viscous_force_N);
end

function dydt = kartODE(~,y,p,useCoulomb)
v = max(y(2),0);
slope = hillSlope(y(1),p);
a = localForces(y(1),v,p,useCoulomb);
% Smooth trial steps through zero speed; terminal event truncates the run.
dydt = [v/sqrt(1+slope^2);a;v];
end

function [a,angleDeg,Fgravity,Frolling,Fdrag,Fviscous] = localForces(x,v,p,useCoulomb)
slope = hillSlope(x,p);
q = sqrt(1+slope.^2);
angleDeg = -atan(slope)*180/pi;
Fgravity = -p.mass*9.81*slope./q;
[Frolling,Fviscous] = frictionForces(v,q,p,useCoulomb);
Fdrag = p.dragFactor*v.^2;
% Retain the uploaded dynamic model's acceleration scaling.
a = ((Fgravity-Frolling-Fdrag-Fviscous)/p.mEffective);
end

function [Frolling,Fviscous] = frictionForces(v,q,p,useCoulomb)
coefficient = p.Crr.*useCoulomb + p.Crr2.*~useCoulomb;
Frolling = coefficient.*p.mass*9.81./q;
speed = max(v,0);
Fviscous = p.cViscous1.*speed.*useCoulomb ...
         + p.cViscous.*speed.*(~useCoulomb);
end

function [value,isterminal,direction] = stopEvents(t,y,p,xTarget)
stopValue = y(2);
% Starting from rest with positive acceleration is handled separately from
% a later downward zero-speed crossing. Avoid a spurious root at release.
if t == 0 && p.v0 == 0, stopValue = 1; end
value = [xTarget-y(1);stopValue];
isterminal = [1;1];
direction = [-1;-1];
end

function values = struct2arrayLocal(s)
names = fieldnames(s);
values = zeros(1,numel(names));
for k = 1:numel(names), values(k) = s.(names{k}); end
end


%% Independent work-energy method (distance is the independent variable)
function [results,status] = simulateEnergy(p)
% K includes translation and wheel rotation. With q = sqrt(1+h'(x)^2):
%   dK/dx = -m*g*h'(x) - (FCoulomb + Fviscous + Fdrag)*q
% Energy state = [kinetic energy; aerodynamic work loss; viscous work loss].
% Fviscous = b*v, dWviscous/dx = b*v*q; power = b*v^2.
% Coulomb work simplifies because N~m*g/q and ds=q*dx:
%   WCoulomb = m*g*(mu1*(min(x,xFrictionSwitch)-xStart)
%                      + mu2*max(x-xFrictionSwitch,0)).
% This solver does not use dynamic speeds, times, or accelerations.
K0 = 0.5*p.mEffective*p.v0^2;
initialSlope = hillSlope(p.xStart,p);
initialGradient = -p.mass*9.81*(initialSlope+p.Crr);
if K0 == 0 && initialGradient <= 0
    x_m = p.xStart;
    state = [0 0 0];
    status = 'Energy method: remained at rest.';
else
    options = odeset('RelTol',1e-9,'AbsTol',[1e-7 1e-7 1e-7], ...
        'MaxStep',min(0.05,(p.xEnd-p.xStart)/1000), ...
        'Events',@(x,z) energyStop(x,z,p,K0));
    [x_m,state,~,~,eventIndex] = ode45(@(x,z) energyODE(x,z,p,true), ...
        [p.xStart p.xFrictionSwitch],[K0;0;0],options);
    if isempty(eventIndex)
        [viscousX,viscousState,~,~,eventIndex] = ode45(@(x,z) energyODE(x,z,p,false), ...
            [p.xFrictionSwitch p.xEnd],state(end,:).',options);
        x_m = [x_m;viscousX(2:end)];
        state = [state;viscousState(2:end,:)];
    end
    if isempty(eventIndex)
        status = 'Energy method: reached end x.';
    else
        state(end,1) = 0;
        status = 'Energy method: stopped where kinetic energy reached zero.';
    end
end
kinetic_J = max(state(:,1),0);
speed_mps = sqrt(2*kinetic_J/p.mEffective);
elevation_m = hillHeight(x_m,p);
original_x_m = x_m+p.sourceXOffset; % Source coordinate of the hill polynomial.
horizontal_from_start_m = x_m-p.xStart;
original_elevation_m = elevation_m+p.sourceElevationAtStart;
initial_kinetic_J = repmat(K0,numel(x_m),1);
gravity_work_J = p.mass*9.81*(hillHeight(p.xStart,p)-elevation_m);
rolling_loss_J = p.mass*9.81*( ...
    p.Crr*(min(x_m,p.xFrictionSwitch)-p.xStart) + ...
    p.Crr2*max(x_m-p.xFrictionSwitch,0));
drag_loss_J = state(:,2);
viscous_loss_J = state(:,3);
% Conservation diagnostic: should remain near zero, including with losses.
balance_residual_J = kinetic_J-initial_kinetic_J-gravity_work_J+rolling_loss_J+drag_loss_J+viscous_loss_J;
results = table(x_m,original_x_m,horizontal_from_start_m,elevation_m,original_elevation_m,speed_mps, ...
    kinetic_J,initial_kinetic_J,gravity_work_J,rolling_loss_J,drag_loss_J,viscous_loss_J,balance_residual_J);
end

function dzdx = energyODE(x,z,p,useCoulomb)
slope = hillSlope(x,p);
q = sqrt(1+slope^2);
speedSquared = 2*max(z(1),0)/p.mEffective;
dragWorkGradient = p.dragFactor*speedSquared*q;
[Frolling,Fviscous] = frictionForces(sqrt(speedSquared),q,p,useCoulomb);
viscousWorkGradient = Fviscous*q;
dzdx = [-p.mass*9.81*slope-Frolling*q-dragWorkGradient-viscousWorkGradient; ...
    dragWorkGradient;viscousWorkGradient];
end

function [value,isterminal,direction] = energyStop(x,z,p,K0)
value = z(1);
if x == p.xStart && K0 == 0, value = 1; end
isterminal = 1;
direction = -1;
end

function sampled = sampleAtX(x,values,query)
% Interpolate kinetic energy before taking sqrt, so release from rest is
% represented accurately. Never extrapolate past an energy-method stop.
sampled = nan(size(query));
tolerance = 1e-7;
if numel(x) == 1
    sampled(abs(query-x(1)) <= tolerance) = values(1);
else
    valid = query >= x(1)-tolerance & query <= x(end)+tolerance;
    clipped = min(max(query(valid),x(1)),x(end));
    sampled(valid) = interp1(x,values,clipped,'pchip');
end
end


function slope = hillSlope(x,p)
% Reduce spatial downhill-angle variation exactly 50%, keeping its start value.
originalAngle = -atan(polyval(p.slopeCoefficients,x));
angle = p.startAngle_rad+p.angleVariationScale*(originalAngle-p.startAngle_rad);
slope = -tan(angle);
end

function elevation = hillHeight(x,p)
% Integrate exactly the same slope used for gravity and path geometry.
% Accumulate between sorted query points to avoid repeated long integrals.
shape = size(x);
[nodes,~,mapping] = unique([p.xStart;x(:)]);
heights = zeros(size(nodes));
anchor = find(nodes==p.xStart,1);
for k=anchor+1:numel(nodes)
    heights(k)=heights(k-1)+integral(@(u) hillSlope(u,p),nodes(k-1),nodes(k), ...
        'AbsTol',1e-10,'RelTol',1e-10);
end
for k=anchor-1:-1:1
    heights(k)=heights(k+1)-integral(@(u) hillSlope(u,p),nodes(k),nodes(k+1), ...
        'AbsTol',1e-10,'RelTol',1e-10);
end
elevation = reshape(heights(mapping(2:end)),shape);
end

function shifted = shiftPolynomial(coefficients,offset)
% Horner composition gives coefficients of P(x+offset), without refitting.
shifted = coefficients(1);
for k = 2:numel(coefficients)
    shifted = conv(shifted,[1 offset]);
    shifted(end) = shifted(end)+coefficients(k);
end
end


function [xVisible,yVisible] = trimComparisonSamples(x,y)
% Clip connected comparison curves at x=0, interpolating the boundary.
% Input x is nondecreasing; y may contain speed and both band edges.
firstVisible = find(x >= 0,1,'first');
if isempty(firstVisible)
    % Keep a graphics handle for checkboxes, even when the curve is empty.
    % Selected CSV export removes these placeholders.
    xVisible = NaN;
    yVisible = nan(1,size(y,2));
    return
end
xVisible = x(firstVisible:end);
yVisible = y(firstVisible:end,:);
if firstVisible > 1 && x(firstVisible) > 0
    before = firstVisible-1;
    fraction = -x(before)/(x(firstVisible)-x(before));
    yAtZero = y(before,:)+fraction*(y(firstVisible,:)-y(before,:));
    xVisible = [0;xVisible];
    yVisible = [yAtZero;yVisible];
end
end


%% Embedded averaged GPS speed (no external files required)
function gps = embeddedGPSSpeed()
% Source: AveragedKartSpeed_0_to_30m.csv, columns x_m and speed_mps.
% Exact supplied samples; no smoothing, shifts, or speed fit applied.
% The CSV run1/run2 and lower/upper envelope columns are not plotted.
% The selectable GPS band uses the existing assumed +/-15% speed error.
% Columns: comparison coordinate (m), averaged GPS speed (m/s).
data = [
    0.0 0.0;
    1.0 0.3931;
    2.0 0.6236;
    3.0 0.8834;
    4.0 1.1183;
    5.0 1.295;
    6.0 1.4325;
    7.0 1.5654;
    8.0 1.6838;
    9.0 1.7688;
    10.0 1.8471;
    11.0 1.9383;
    12.0 2.0287;
    13.0 2.1146;
    14.0 2.2006;
    15.0 2.2972;
    16.0 2.3932;
    17.0 2.488;
    18.0 2.5728;
    19.0 2.6588;
    20.0 2.7456;
    21.0 2.8285;
    22.0 2.8996;
    23.0 2.9707;
    24.0 3.0427;
    25.0 3.1173;
    26.0 3.1924;
    27.0 3.2675;
    28.0 3.3284;
    29.0 3.39;
    30.0 3.4516;
    ];
gps = table(data(:,1),data(:,2),'VariableNames',{'x_m','speed_mps'});
end
