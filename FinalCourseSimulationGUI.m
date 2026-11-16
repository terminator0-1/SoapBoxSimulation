function FinalCourseSimulationLagrange


spec = {
    'mass', 'Kart + rider mass (kg)', 117, [0.001 Inf], 'Total mass INCLUDING all wheels.';
    'v0', 'Initial velocity (m/s)', 0, [0 Inf], 'Speed along the hill at x = 0 m, used by all three solvers.';
    'tMax', 'Time limit (s)', 120, [0.01 3600], 'Time limit for dynamic and Lagrangian solvers.';
    'pushForce_N', 'Start push over 1 m (N)', 200, [0 Inf], 'Constant forward force along the first 1 m of travel. Set 0 to disable.';
    'Crr1', 'Coulomb coefficient, μ₁', 0.02, [0 Inf], 'Section 1: 0 <= x < 31 m. Surface obstacles override this locally.';
    'b1', 'Linear b₁ (N*s/m)', 18, [0 Inf], 'Section 1: Fviscous = b1*v. Surface obstacles do not change b.';
    'Crr2', 'Coulomb coefficient, μ₂', 0.02, [0 Inf], 'Section 2: 31 <= x < 93 m. Surface obstacles override this locally.';
    'b2', 'Linear b₂ (N*s/m)', 18, [0 Inf], 'Section 2: Fviscous = b2*v. Surface obstacles do not change b.';
    'Crr3', 'Coulomb coefficient, μ₃', 0.02, [0 Inf], 'Section 3: x >= 93 m. Surface obstacles override this locally.';
    'b3', 'Linear b₃ (N*s/m)', 18, [0 Inf], 'Section 3: Fviscous = b3*v. Surface obstacles do not change b.'
    };

fixed = struct('xStart',0,'xEnd',116,'pushDistance_m',1,'sourceXOffset',0,'sectionBoundary1_m',31,'sectionBoundary2_m',93, ...
    'CdA',0,'rho',1.225,'nWheels',4,'radius',0.127,'Iwheel',0.0175, ...
    'speedErrorFraction',0.15,'simulatedSpeedOffset_mps',0, ...
    'comparisonXOffset_m',0.4,'angleSmoothingFraction',0,'lagrangianAngleSmoothingFraction',0.30,'angleSmoothingSigma_m',10.5);
nInputs = size(spec,1);
nRiderInputs = 4;
fields = struct();
obstacles = emptyObstacles();
lastResults = table();
lastParameters = struct();
lastLagrangianParameters = struct();
lastEnergy = table();
lastLagrangian = table();
gpsSpeed = embeddedGPSSpeed(fixed.xEnd);
busy = false;
ink = [0.12 0.19 0.27];
blue = [0.08 0.43 0.72];
green = [0.10 0.58 0.43];
orange = [0.90 0.43 0.12];
animationElapsed = 0;
animationBase = 0;
animationClock = [];
animationTimer = [];
animationPlots = struct();


fig = uifigure('Name','Downhill Kart | Dynamics, Energy & Lagrangian', ...
    'Position',[45 35 1280 870], 'Color',[0.96 0.97 0.98], 'Visible','off');
outer = uigridlayout(fig,[4 2]);
outer.RowHeight = {40,'1x',72,30};
outer.ColumnWidth = {365,'1x'};
outer.Padding = [12 8 12 8];
outer.ColumnSpacing = 12;
outer.RowSpacing = 7;
header = uigridlayout(outer,[1 3]);
header.Layout.Row = 1; header.Layout.Column = [1 2];
header.ColumnWidth = {'1x',145,145};
header.Padding = [0 0 0 0];
uilabel(header,'Text','DOWNHILL / Three-model comparison', ...
    'FontSize',23,'FontWeight','bold','FontColor',ink);
uibutton(header,'Text','Import settings','ButtonPushedFcn',@importSettings);
uibutton(header,'Text','Export settings','ButtonPushedFcn',@exportSettings);
leftTabs = uitabgroup(outer);
leftTabs.Layout.Row = 2; leftTabs.Layout.Column = 1;
panel = uitab(leftTabs,'Title','Rider & start');
controls = uigridlayout(panel,[nRiderInputs+7 2]);
controls.ColumnWidth = {'1x',100};
controls.RowHeight = [repmat({38},1,nRiderInputs), {28,38,32,32,138,32,'1x'}];
controls.Padding = [12 8 12 8];
controls.RowSpacing = 5;
controls.Scrollable = 'on';
frictionTab = uitab(leftTabs,'Title','Friction sections');
frictionGrid = uigridlayout(frictionTab,[4 1]);
frictionGrid.RowHeight = {118,118,118,'1x'};
frictionGrid.Padding = [10 10 10 10];
frictionGrid.Scrollable = 'on';
sectionTitles = {'Section 1 | 0-31 m','Section 2 | 31-93 m','Section 3 | 93 m to finish'};
sectionGrids = cell(1,3);
for section = 1:3
    sectionPanel = uipanel(frictionGrid,'Title',sectionTitles{section},'FontWeight','bold');
    sectionGrids{section} = uigridlayout(sectionPanel,[2 2]);
    sectionGrids{section}.ColumnWidth = {'1x',95};
    sectionGrids{section}.RowHeight = {34,34};
    sectionGrids{section}.Padding = [8 6 8 6];
end
uilabel(frictionGrid,'Text',sprintf([ ...
    'F = mu(x)*m*g*cos(theta) + b(x)*v\n\n' ...
    'Sections use hill distance, before comparison-plot offsets.\n\n' ...
    'Surface obstacles replace mu only. The section b still applies.\n\n' ...
    'At x = 31 m section 2 starts; at x = 93 m section 3 starts.']), ...
    'WordWrap','on','VerticalAlignment','top','FontSize',12);
for i = 1:size(spec,1)
    if i <= nRiderInputs
        inputParent = controls; inputRow = i;
    else
        section = ceil((i-nRiderInputs)/2);
        inputParent = sectionGrids{section};
        inputRow = mod(i-nRiderInputs-1,2)+1;
    end
    label = uilabel(inputParent,'Text',spec{i,2},'FontSize',12);
    label.Layout.Row = inputRow; label.Layout.Column = 1;
    f = uieditfield(inputParent,'numeric','Value',spec{i,3}, ...
        'Limits',spec{i,4},'Tooltip',spec{i,5}, ...
        'ValueChangedFcn',@parametersChanged);
    f.Layout.Row = inputRow; f.Layout.Column = 2;
    fields.(spec{i,1}) = f;
end

autoRun = uicheckbox(controls,'Text','Update plots after each edit', ...
    'Value',true,'ValueChangedFcn',@parametersChanged);
autoRun.Layout.Row = nRiderInputs+1; autoRun.Layout.Column = [1 2];
runButton = uibutton(controls,'Text','Run simulation','FontWeight','bold', ...
    'BackgroundColor',blue,'FontColor',[1 1 1], 'ButtonPushedFcn',@runSimulation);
runButton.Layout.Row = nRiderInputs+2; runButton.Layout.Column = [1 2];
resetButton = uibutton(controls,'Text','Reset all parameters','ButtonPushedFcn',@resetInputs);
resetButton.Layout.Row = nRiderInputs+3; resetButton.Layout.Column = [1 2];
exportButton = uibutton(controls,'Text','Export results to CSV', ...
    'Enable','off','ButtonPushedFcn',@exportResults);
exportButton.Layout.Row = nRiderInputs+4; exportButton.Layout.Column = [1 2];
note = uilabel(controls,'Text',sprintf(['Simulation: x = 0-116 m; extra 30%% smoothing only for Lagrangian.\n' ...
    'Resistance: F = b(x)*v + mu(x)*m*g*cos(theta).\n' ...
    'Edit mu and b in the Friction sections tab. Surfaces override only mu.\n' ...
    'Push acts over the first 1 m of travel (0 N disables).\n' ...
    'Aerodynamic drag is neglected.']), ...
    'WordWrap','on','FontColor',[0.35 0.39 0.44],'VerticalAlignment','top');
note.Layout.Row = nRiderInputs+5; note.Layout.Column = [1 2];
energyExportButton = uibutton(controls,'Text','Export full energy CSV', ...
    'Enable','off','ButtonPushedFcn',@exportEnergy);
energyExportButton.Layout.Row = nRiderInputs+6; energyExportButton.Layout.Column = [1 2];

obstacleTab = uitab(leftTabs,'Title','Course obstacles');


    editor = uigridlayout(obstacleTab, [17 2]);
    editor.ColumnWidth = {'1x', 130};
    editor.Padding = [12 8 12 8];
    editor.RowSpacing = 5;
    editor.RowHeight = {22, 100, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 58, 32, 32, 40, 85};
    editor.Scrollable = 'on';
    courseLengthLabel = addLabel(editor, 1, sprintf('Course length: %.2f m', fixed.xEnd), true);
    obstacleList = uilistbox(editor, 'Items', {'[New obstacle]'}, ...
                             'Value', '[New obstacle]', 'ValueChangedFcn', @(~, ~) selectObstacle());
    place(obstacleList, 2, [1 2]);
    l = addLabel(editor, 3, 'Obstacle type', false);
    l.Layout.Column = 1;
    typeField = uidropdown(editor, 'Items', {'Ramp', 'Surface'}, 'Value', 'Ramp', ...
                           'ValueChangedFcn', @(~, ~) updateEditor());
    place(typeField, 3, 2);
    l = addLabel(editor, 4, 'Name (optional)', false);
    l.Layout.Column = 1;
    nameField = uieditfield(editor, 'text', 'Value', '');
    place(nameField, 4, 2);
    startField = numericControl(editor, 5, 'Start station (m)', min(20, fixed.xEnd / 3), [0 fixed.xEnd]);
    heightField = numericControl(editor, 6, 'Ramp height above terrain (m)', 0.3, [0 20]);
    upField = numericControl(editor, 7, 'Uphill horizontal length (m)', min(2, fixed.xEnd / 10), [0 fixed.xEnd]);
    flatField = numericControl(editor, 8, 'Top-section horizontal length (m)', 0, [0 fixed.xEnd]);
    downField = numericControl(editor, 9, 'Downhill horizontal length (m)', min(2, fixed.xEnd / 10), [0 fixed.xEnd]);
    l = addLabel(editor, 10, 'Surface name', false);
    l.Layout.Column = 1;
    surfaceField = uidropdown(editor, 'Items', {'Custom', 'Asphalt', 'Concrete', 'Grass', 'Gravel', 'Sand', 'Wood'}, ...
                              'Value', 'Custom');
    place(surfaceField, 10, 2);
    muField = numericControl(editor, 11, 'Surface coefficient μ', 0.02, [0 1]);
    lengthField = numericControl(editor, 12, 'Surface horizontal length (m)', min(5, fixed.xEnd / 5), [0 fixed.xEnd]);
    editorInfo = addLabel(editor, 13, '', false);
    editorInfo.WordWrap = 'on';
    editorInfo.FontSize = 11;
    addButton = uibutton(editor, 'Text', 'Add as new', 'ButtonPushedFcn', @(~, ~) commitObstacle(false));
    place(addButton, 14, 1);
    applyButton = uibutton(editor, 'Text', 'Apply edits', 'ButtonPushedFcn', @(~, ~) commitObstacle(true));
    place(applyButton, 14, 2);
    deleteButton = uibutton(editor, 'Text', 'Delete selected', 'ButtonPushedFcn', @deleteObstacle);
    place(deleteButton, 15, [1 2]);
    editorMessage = addLabel(editor, 16, 'Select New, fill in parameters, then Add.', false);
    editorMessage.WordWrap = 'on';
    editorMessage.FontSize = 11;
    l = addLabel(editor, 17, ...
                 'Positions start at the measured terrain (station 0). A surface may cover a ramp. Overlapping ramps or overlapping surface sections are not allowed.', false);
    l.WordWrap = 'on';
    l.FontSize = 11;
    for c = {startField, heightField, upField, flatField, downField, lengthField}
        c{1}.ValueChangedFcn = @(~, ~) updateEditor();
    end


metrics = uigridlayout(outer,[1 3]);
metrics.Layout.Row = 3; metrics.Layout.Column = [1 2];
metrics.Padding = [12 6 12 6];
metrics.BackgroundColor = [1 1 1];
metricSpeed = uilabel(metrics,'Text','Final speed','FontSize',15,'FontWeight','bold','FontColor',blue);
metricTime = uilabel(metrics,'Text','Elapsed time','FontSize',15,'FontWeight','bold','FontColor',ink);
metricMass = uilabel(metrics,'Text','Effective mass','FontSize',15,'FontWeight','bold','FontColor',ink);
tabs = uitabgroup(outer);
tabs.Layout.Row = 2; tabs.Layout.Column = 2;
fullTab = uitab(tabs,'Title','Full course');
fullGrid = uigridlayout(fullTab,[4 1]);
fullGrid.RowHeight = {'1x','1x','1x',105};
fullGrid.Padding = [4 2 4 2]; fullGrid.RowSpacing = 5;
axFullHill = uiaxes(fullGrid);
axFullSpeed = uiaxes(fullGrid);
axFullEnergy = uiaxes(fullGrid);
axFullForce = uiaxes(fullGrid);
linkaxes([axFullHill axFullSpeed axFullEnergy axFullForce],'x');

animationTab = uitab(tabs,'Title','Animation');
animationLayout = uigridlayout(animationTab,[3 1]);
animationLayout.RowHeight = {36,30,'1x'};
animationLayout.Padding = [8 8 8 8];
animToolbar = uigridlayout(animationLayout,[1 4]);
animToolbar.ColumnWidth = {85,85,85,'1x'};
animToolbar.Padding = [0 0 0 0];
playButton = uibutton(animToolbar,'Text','Play','Enable','off','ButtonPushedFcn',@playAnimation);
pauseButton = uibutton(animToolbar,'Text','Pause','Enable','off','ButtonPushedFcn',@pauseAnimation);
uibutton(animToolbar,'Text','Restart','ButtonPushedFcn',@restartAnimation);
animationStatus = uilabel(animToolbar,'Text','Ready | real time (1x)');
animationReadout = uilabel(animationLayout,'Text','', ...
    'FontSize',16,'FontWeight','bold','FontColor',blue);
animationPanel = uipanel(animationLayout,'BorderType','none');
animSpeedAxes = uiaxes(animationPanel,'Units','normalized', ...
    'Position',[0.11 0.59 0.85 0.34],'PositionConstraint','innerposition');
animHillAxes = uiaxes(animationPanel,'Units','normalized', ...
    'Position',[0.11 0.11 0.85 0.34],'PositionConstraint','innerposition');
linkaxes([animSpeedAxes animHillAxes],'x');
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
comparisonLabels = {'Dynamic','Energy balance','Real kart speed', ...
    'Dynamic +/-15% band','Real kart +/-15% band','Lagrangian'};
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
    'Simulation: physical speed; comparison x = hill x + 0.4 m.\n' ...
    'RealKartData: left shift 19.07 m; vertical start; clipped at 116 m.']), ...
    'WordWrap','on','FontSize',11);
bandNote.Layout.Row = 1;
bandNote.Layout.Column = 1;
motionTab = uitab(tabs,'Title','Motion');
hillTab = uitab(tabs,'Title','Hill profile');
energyTab = uitab(tabs,'Title','Energy balance');
energyGrid = uigridlayout(energyTab,[4 1]);
energyGrid.RowHeight = {66,'1x','1x',60};
energyNote = uilabel(energyGrid,'Text',sprintf([ ...
    'K = K0 + m*g*(h0-h) + Wpush - WCoulomb - Wviscous;  K = 0.5*mEffective*v^2\n' ...
    'Linear b(x)*v with separate mu and b in 0-31, 31-93, and 93 m onward.\n' ...
    'Course: 0 to 116 m (smoothed merged hill profile). Time limit applies to dynamic and Lagrangian solvers.']), ...
    'WordWrap','on','FontSize',12);
axEnergy = uiaxes(energyGrid);
axAgreement = uiaxes(energyGrid);
energySummary = uilabel(energyGrid,'Text','Ready.','WordWrap','on');
lagrangianTab = uitab(tabs,'Title','Lagrangian');
lagrangianGrid = uigridlayout(lagrangianTab,[5 1]);
lagrangianGrid.RowHeight = {90,'1x','1x',52,34};
uilabel(lagrangianGrid,'Text',sprintf([ ...
    'L = 0.5*mEffective*(1+h''(x)^2)*xdot^2 - m*g*h(x)\n' ...
    'mEffective*((1+h''^2)*xddot + h''*h''''*xdot^2) + m*g*h'' = Qx\n' ...
    'Qx = (Fpush - FCoulomb - b*v - Fdrag)*sqrt(1+h''^2); v = xdot*sqrt(1+h''^2)\n' ...
    'Lagrangian only: extra 30%% angle smoothing. Same wheel inertia, forces, and obstacles.']), ...
    'WordWrap','on','FontSize',12);
axLagrangianSpeed = uiaxes(lagrangianGrid);
axLagrangianResidual = uiaxes(lagrangianGrid);
lagrangianSummary = uilabel(lagrangianGrid,'Text','Ready.','WordWrap','on');
lagrangianExportButton = uibutton(lagrangianGrid,'Text','Export full Lagrangian CSV', ...
    'Enable','off','ButtonPushedFcn',@exportLagrangian);
plotGrid = uigridlayout(motionTab,[2 2]);
plotGrid.Padding = [0 0 0 0];
axSpeed = uiaxes(plotGrid); axSpeed.Layout.Row = 1; axSpeed.Layout.Column = 1;
axDistanceSpeed = uiaxes(plotGrid); axDistanceSpeed.Layout.Row = 1; axDistanceSpeed.Layout.Column = 2;
axAccel = uiaxes(plotGrid); axAccel.Layout.Row = 2; axAccel.Layout.Column = 1;
axDistance = uiaxes(plotGrid); axDistance.Layout.Row = 2; axDistance.Layout.Column = 2;
hillGrid = uigridlayout(hillTab,[3 1]);
hillGrid.RowHeight = {48,'1x','1x'};
formulaLabel = uilabel(hillGrid,'Text', ...
    'Dynamic / energy: original terrain smoothing. Lagrangian: extra 30% angle-smoothing blend. Both retain the total height drop.', ...
    'WordWrap','on','FontSize',12);
axHill = uiaxes(hillGrid);
axGrade = uiaxes(hillGrid);
statusLabel = uilabel(outer,'Text','Ready.','WordWrap','on', ...
    'FontColor',[0.30 0.34 0.40]);
statusLabel.Layout.Row = 4; statusLabel.Layout.Column = [1 2];
statusLabel.FontSize = 11;
for ax = [axFullHill axFullSpeed axFullEnergy axFullForce axComparison ...
        axSpeed axDistanceSpeed axAccel axDistance axHill axGrade axEnergy ...
        axAgreement animSpeedAxes animHillAxes]
    ax.FontSize = 10;
    ax.XColor = ink; ax.YColor = ink;
    ax.GridAlpha = 0.12;
    grid(ax,'on'); box(ax,'off');
end
animationTimer = timer('ExecutionMode','fixedRate','Period',0.033, ...
    'BusyMode','drop','Name','Kart dynamics playback', ...
    'TimerFcn',@animationTick,'ErrorFcn',@animationError);
fig.DeleteFcn = @cleanupAnimation;
tabs.SelectionChangedFcn = @tabChanged;
refreshList(0);
runSimulation();
fig.Visible = 'on';


    function index = selectedIndex()
        index = find(strcmp(obstacleList.Items, obstacleList.Value), 1) - 1;
        if isempty(index)
            index = 0;
        end
    end

    function updateEditor()
        isRamp = strcmp(typeField.Value, 'Ramp');
        for field = {heightField, upField, flatField, downField}
            field{1}.Enable = onOff(isRamp);
        end
        for field = {surfaceField, muField, lengthField}
            field{1}.Enable = onOff(~isRamp);
        end
        applyButton.Enable = onOff(selectedIndex() > 0);
        deleteButton.Enable = onOff(selectedIndex() > 0);
        if isRamp
            finish = startField.Value + upField.Value + flatField.Value + downField.Value;
            editorInfo.Text = sprintf('Ends at %.2f m. Height is added above terrain; the top follows terrain grade. No airborne jump.', finish);
        else
            editorInfo.Text = sprintf('Ends at %.2f m. This section overrides mu only; the local friction-section b still applies.', startField.Value + lengthField.Value);
        end
    end

    function selectObstacle()
        index = selectedIndex();
        if index == 0
            o = newObstacle();
            o.start_m = min(20, fixed.xEnd / 3);
            o.up_m = min(2, fixed.xEnd / 10);
            o.down_m = o.up_m;
            o.length_m = min(5, fixed.xEnd / 5);
            o.mu = min(fields.Crr1.Value,1);
        else
            o = obstacles(index);
        end
        setEditorFields(o);
        editorMessage.Text = 'Change fields, then click Add as new or Apply edits.';
    end

    function o = getEditorFields()
        o = newObstacle();
        o.type = typeField.Value;
        o.name = nameField.Value;
        o.start_m = startField.Value;
        o.height_m = heightField.Value;
        o.up_m = upField.Value;
        o.flat_m = flatField.Value;
        o.down_m = downField.Value;
        o.surface = surfaceField.Value;
        o.mu = muField.Value;
        o.length_m = lengthField.Value;
    end

    function setEditorFields(o)
        typeField.Value = o.type;
        nameField.Value = o.name;
        startField.Value = o.start_m;
        heightField.Value = o.height_m;
        upField.Value = o.up_m;
        flatField.Value = o.flat_m;
        downField.Value = o.down_m;
        surfaceField.Value = o.surface;
        muField.Value = o.mu;
        lengthField.Value = o.length_m;
        updateEditor();
    end

    function refreshList(index)
        items = cell(1, numel(obstacles) + 1);
        items{1} = '[New obstacle]';
        for j = 1:numel(obstacles)
            o = obstacles(j);
            items{j + 1} = sprintf('%d | %s | %.1f-%.1f m', j, o.name, o.start_m, obstacleEnd(o));
        end
        obstacleList.Items = items;
        obstacleList.Value = items{index + 1};
        selectObstacle();
    end

    function commitObstacle(editExisting)
        index = selectedIndex();
        if editExisting && index == 0
            return
        end
        o = getEditorFields();
        o.name = strtrim(o.name);
        candidate = obstacles;
        if ~editExisting
            index = numel(candidate) + 1;
        end
        if isempty(o.name)
            o.name = sprintf('%s %d', o.type, index);
        end
        candidate(index) = o;
        try
            validateObstacles(candidate, fixed.xEnd);
        catch err
            uialert(fig, err.message, 'Check obstacle parameters');
            return
        end
        obstacles = candidate;
        refreshList(index);
        runSimulation();
        editorMessage.Text = sprintf('Saved %s at %.2f-%.2f m.', o.name, o.start_m, obstacleEnd(o));
    end

    function deleteObstacle(~, ~)
        index = selectedIndex();
        if index == 0
            return
        end
        obstacles(index) = [];
        refreshList(0);
        runSimulation();
        editorMessage.Text = 'Obstacle removed.';
    end


    function exportSettings(~,~)
        [name,folder] = uiputfile('*.mat','Export settings','KartDynamicsSettings.mat');
        if isequal(name,0), return; end
        settings = struct();
        for j = 1:nInputs
            settings.(spec{j,1}) = fields.(spec{j,1}).Value;
        end
        settings.obstacles = obstacles;
        settings.editorDraft = getEditorFields();
        settings.selectedObstacle = selectedIndex();
        try
            save(fullfile(folder,name),'settings');
            statusLabel.Text = 'Rider and resistance settings exported.';
        catch err
            uialert(fig,err.message,'Export settings');
        end
    end

    function importSettings(~,~)
        [name,folder] = uigetfile('*.mat','Import kart dynamics settings');
        if isequal(name,0), return; end
        try
            loaded = load(fullfile(folder,name),'settings');
            if ~isfield(loaded,'settings') || ~isstruct(loaded.settings) || ~isscalar(loaded.settings)
                error('Choose a settings file exported by this dynamics app.');
            end


            for section = 1:3
                muKey = sprintf('Crr%d',section);
                bKey = sprintf('b%d',section);
                if ~isfield(loaded.settings,muKey) && isfield(loaded.settings,'Crr')
                    loaded.settings.(muKey) = loaded.settings.Crr;
                end
                if ~isfield(loaded.settings,bKey) && isfield(loaded.settings,'cViscous')
                    loaded.settings.(bKey) = loaded.settings.cViscous;
                end
            end
            if ~isfield(loaded.settings,'pushForce_N')
                loaded.settings.pushForce_N = 0;
            end
            values = zeros(nInputs,1);
            for j = 1:nInputs
                key = spec{j,1};
                if ~isfield(loaded.settings,key), error('Missing setting: %s.',key); end
                value = loaded.settings.(key);
                limits = spec{j,4};
                if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
                        ~isfinite(value) || value < limits(1) || value > limits(2)
                    error('Invalid value for %s.',key);
                end
                values(j) = value;
            end
            candidate = emptyObstacles();
            if isfield(loaded.settings,'obstacles'), candidate = loaded.settings.obstacles; end
            validateObstacles(candidate,fixed.xEnd);
            for k=1:numel(candidate), validateEditorDraft(candidate(k),fixed.xEnd); end
            candidate = orderfields(candidate,newObstacle());
            draft = [];
            if isfield(loaded.settings,'editorDraft')
                draft = loaded.settings.editorDraft;
                validateEditorDraft(draft,fixed.xEnd);
            end
            selected = 0;
            if isfield(loaded.settings,'selectedObstacle')
                selected = loaded.settings.selectedObstacle;
                validateattributes(selected,{'numeric'},{'scalar','real','finite','integer','>=',0,'<=',numel(candidate)});
            end
            for j = 1:nInputs, fields.(spec{j,1}).Value = values(j); end
            obstacles = candidate;
            refreshList(selected);
            if ~isempty(draft), setEditorFields(draft); end
            runSimulation();
        catch err
            uialert(fig,err.message,'Import settings');
        end
    end

    function updateCourseOverview(r,e,l,p)
        xx = unique([linspace(p.xStart,p.xEnd,800),obstacleBreaks(p)]);
        hh = hillHeight(xx,p)+p.sourceElevationAtStart;
        plotResult(axFullHill,xx,hh,'Merged hill profile','','Elevation (m)');
        hold(axFullHill,'on');
        plot(axFullHill,r.x_m,r.original_elevation_m,'Color',orange,'LineWidth',2);
        plot(axFullHill,r.x_m(end),r.original_elevation_m(end),'o', ...
            'MarkerFaceColor',orange,'MarkerEdgeColor',ink,'MarkerSize',7);
        hold(axFullHill,'off');
        plotResult(axFullSpeed,r.x_m,r.speed_mps,'Velocity along the course','','Speed (m/s)');
        hold(axFullSpeed,'on');
        plot(axFullSpeed,e.x_m,e.speed_mps,'--','Color',green,'LineWidth',1.6);
        plot(axFullSpeed,l.x_m,l.speed_mps,':','Color',[0.58 0.20 0.75],'LineWidth',2.1);
        hold(axFullSpeed,'off');
        legend(axFullSpeed,{'Dynamic','Energy balance','Lagrangian'},'Location','northeast');
        markObstacles(axFullSpeed,p,false);
        plotResult(axFullEnergy,e.x_m,e.kinetic_J,'Energy accounting','','Energy (J)');
        hold(axFullEnergy,'on');
        plot(axFullEnergy,e.x_m,e.initial_kinetic_J+e.gravity_work_J+e.push_work_J, ...
            ':','Color',[.53 .30 .72],'LineWidth',1.5);
        plot(axFullEnergy,e.x_m,e.rolling_loss_J+e.viscous_loss_J, ...
            '--','Color',orange,'LineWidth',1.5);
        hold(axFullEnergy,'off');
        legend(axFullEnergy,{'Kinetic','Initial + gravity + push','Resistance losses'},'Location','northeast');
        plotResult(axFullForce,r.x_m,r.rolling_force_N+r.viscous_force_N, ...
            'Rolling + speed-dependent resistance','Horizontal station from terrain start (m)','Force (N)');
        xlim(axFullHill,[p.xStart p.xEnd]);
        markObstacles(axFullHill,p,true);
        markObstacles(axFullForce,p,false);
    end

    function prepareAnimation()
        pauseAnimation();
        animationElapsed = 0;
        r = lastResults; p = lastParameters;
        xx = unique([linspace(p.xStart,p.xEnd,800),obstacleBreaks(p)]);
        hh = hillHeight(xx,p)+p.sourceElevationAtStart;
        cla(animSpeedAxes); cla(animHillAxes);
        hold(animSpeedAxes,'on'); hold(animHillAxes,'on');
        plot(animSpeedAxes,r.x_m,r.speed_mps,'Color',[.72 .80 .87],'LineWidth',1.5);
        animationPlots.trail = plot(animSpeedAxes,NaN,NaN,'Color',blue,'LineWidth',2.2);
        plot(animHillAxes,xx,hh,'Color',ink,'LineWidth',2);
        animationPlots.speedCursor = xline(animSpeedAxes,r.x_m(1),':','Color',orange);
        animationPlots.hillCursor = xline(animHillAxes,r.x_m(1),':','Color',orange);
        animationPlots.speedPoint = plot(animSpeedAxes,NaN,NaN,'o', ...
            'MarkerSize',9,'MarkerFaceColor',blue,'Color',blue);
        animationPlots.rider = plot(animHillAxes,NaN,NaN,'o', ...
            'MarkerSize',11,'MarkerFaceColor',orange,'Color',ink);
        hold(animSpeedAxes,'off'); hold(animHillAxes,'off');
        title(animSpeedAxes,'Velocity along the course'); ylabel(animSpeedAxes,'Speed (m/s)');
        title(animHillAxes,'Rider position'); ylabel(animHillAxes,'Elevation (m)');
        xlabel(animHillAxes,'Horizontal station from terrain start (m)');
        markObstacles(animHillAxes,p,true);
        xlim(animSpeedAxes,[p.xStart p.xEnd]);
        ylim(animSpeedAxes,[0 max(1,max(r.speed_mps)*1.12)]);
        ylim(animHillAxes,[min(hh)-.5 max(hh)+.5]);
        for ax = [animSpeedAxes animHillAxes]
            ax.Interactions = []; ax.Toolbar.Visible = 'off';
        end
        animSpeedAxes.InnerPosition = [0.11 0.59 0.85 0.34];
        animHillAxes.InnerPosition = [0.11 0.11 0.85 0.34];
        updateAnimationFrame();
        playButton.Enable = 'on';
        animationStatus.Text = 'Ready | real time (1x)';
    end

    function playAnimation(varargin)
        if busy || isempty(lastResults) || isempty(fieldnames(animationPlots)), return; end
        if strcmp(animationTimer.Running,'on'), return; end
        if animationElapsed >= lastResults.time_s(end), animationElapsed = 0; end
        animationBase = animationElapsed;
        animationClock = tic;
        animationStatus.Text = 'Playing | real time (1x)';
        playButton.Enable = 'off'; pauseButton.Enable = 'on';
        start(animationTimer);
    end

    function pauseAnimation(varargin)
        if ~isempty(animationTimer) && isvalid(animationTimer) && strcmp(animationTimer.Running,'on')
            animationElapsed = min(lastResults.time_s(end),animationBase+toc(animationClock));
            stop(animationTimer);
            updateAnimationFrame();
            animationStatus.Text = 'Paused | real time (1x)';
        end
        pauseButton.Enable = 'off';
        if ~isempty(lastResults), playButton.Enable = 'on'; end
    end

    function restartAnimation(varargin)
        pauseAnimation();
        animationElapsed = 0;
        if isempty(lastResults) || isempty(fieldnames(animationPlots)), return; end
        updateAnimationFrame();
        animationStatus.Text = 'Ready | real time (1x)';
    end

    function animationTick(varargin)
        if ~isvalid(fig), return; end
        animationElapsed = min(lastResults.time_s(end),animationBase+toc(animationClock));
        updateAnimationFrame();
        if animationElapsed >= lastResults.time_s(end)
            stop(animationTimer);
            playButton.Enable = 'on'; pauseButton.Enable = 'off';
            animationStatus.Text = 'Finished | real time (1x)';
        end
    end

    function updateAnimationFrame()
        r = lastResults;
        if isempty(r) || isempty(fieldnames(animationPlots)), return; end
        xx = interp1(r.time_s,r.x_m,animationElapsed,'linear');
        vv = interp1(r.time_s,r.speed_mps,animationElapsed,'linear');
        hh = hillHeight(xx,lastParameters)+lastParameters.sourceElevationAtStart;
        mask = r.time_s <= animationElapsed;
        animationPlots.trail.XData = [r.x_m(mask);xx];
        animationPlots.trail.YData = [r.speed_mps(mask);vv];
        animationPlots.speedPoint.XData = xx; animationPlots.speedPoint.YData = vv;
        animationPlots.rider.XData = xx; animationPlots.rider.YData = hh;
        animationPlots.speedCursor.Value = xx; animationPlots.hillCursor.Value = xx;
        animationReadout.Text = sprintf('%.2f m/s  |  %.2f mph  |  x = %.2f m  |  t = %.2f / %.2f s', ...
            vv,vv*2.236936292,xx,animationElapsed,r.time_s(end));
        drawnow limitrate nocallbacks
    end

    function tabChanged(~,~)
        if tabs.SelectedTab ~= animationTab, pauseAnimation(); end
    end

    function animationError(varargin)
        pauseAnimation();
        animationStatus.Text = 'Playback stopped. Press Restart to reset.';
    end

    function cleanupAnimation(varargin)
        if ~isempty(animationTimer) && isvalid(animationTimer)
            stop(animationTimer); delete(animationTimer);
        end
    end

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
        lagrangianExportButton.Enable = 'off';
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
        pauseAnimation();
        playButton.Enable = 'off';
        busy = true;
        runButton.Enable = 'off';
        exportButton.Enable = 'off';
        energyExportButton.Enable = 'off';
        lagrangianExportButton.Enable = 'off';
        comparisonExportButton.Enable = 'off';
        cleanup = onCleanup(@finishRun);
        updateEditor();
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
            validateObstacles(obstacles,p.xEnd);
            p.obstacles = obstacles;
            [result,p,message] = simulateKart(p);
            [energy,energyMessage] = simulateEnergy(p);


            pLag = p;
            pLag.angleSmoothingFraction = p.lagrangianAngleSmoothingFraction;
            pLag = rmfield(pLag,'pushEndX_m');
            [lagHeightPP,~] = mergedHillInterpolants(pLag);
            pLag.sourceElevationAtStart = ppval(lagHeightPP,pLag.sourceXStart)+rampHeight(pLag.xStart,pLag);
            pLag.pushEndX_m = findPushEnd(pLag);
            [lagrangian,lagrangianMessage] = simulateLagrangian(pLag);

            result.energy_kinetic_J = sampleAtX(energy.x_m,energy.kinetic_J,result.x_m);
            result.energy_speed_mps = sqrt(2*max(result.energy_kinetic_J,0)/p.mEffective);
            result.energy_speed_mps(~isfinite(result.energy_kinetic_J)) = NaN;
            result.dynamic_minus_energy_mps = result.speed_mps-result.energy_speed_mps;

            [lagX,lagIndices] = unique(lagrangian.x_m,'stable');
            lagK = sampleAtX(lagX,lagrangian.kinetic_J(lagIndices),result.x_m);
            result.lagrangian_speed_mps = sqrt(2*max(lagK,0)/p.mEffective);
            result.lagrangian_speed_mps(~isfinite(lagK)) = NaN;
            result.lagrangian_minus_dynamic_mps = result.lagrangian_speed_mps-result.speed_mps;
            lagrangian.comparison_x_m = lagrangian.x_m+p.comparisonXOffset_m;
            lagrangian.displayed_speed_mps = lagrangian.speed_mps+p.simulatedSpeedOffset_mps;
            lagrangian.gps_speed_mps = sampleRecordedSpeed(gpsSpeed,lagrangian.comparison_x_m);
            lagrangian.displayed_speed_error_mps = lagrangian.displayed_speed_mps-lagrangian.gps_speed_mps;
            result.lagrangian_displayed_speed_mps = result.lagrangian_speed_mps+p.simulatedSpeedOffset_mps;

            result.comparison_x_m = result.x_m+p.comparisonXOffset_m;
            energy.comparison_x_m = energy.x_m+p.comparisonXOffset_m;


            result.gps_speed_mps = sampleRecordedSpeed(gpsSpeed,result.comparison_x_m);
            energy.gps_speed_mps = sampleRecordedSpeed(gpsSpeed,energy.comparison_x_m);
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
            lastLagrangian = lagrangian;
            lastResults = result;
            lastParameters = p;
            lastLagrangianParameters = pLag;
            fig.UserData = struct('Results',result,'Energy',energy,'Lagrangian',lagrangian,'Parameters',p,'LagrangianParameters',pLag,'GPSSpeed',gpsSpeed);
            plotResult(axSpeed,result.time_s,result.displayed_speed_mps, ...
                'Dynamic velocity over time','Time (s)','Speed (m/s)');
            drawComparison(axDistanceSpeed,result,energy,lagrangian,p,false);
            drawComparison(axComparison,result,energy,lagrangian,p,true);
            plotResult(axAccel,result.time_s,result.acceleration_mps2, ...
                'Acceleration over time','Time (s)','Acceleration (m/s^2)');
            plotResult(axDistance,result.time_s,result.distance_m, ...
                'Distance over time','Time (s)','Distance along hill (m)');
            plotResult(axEnergy,energy.x_m,energy.kinetic_J, ...
                'Work-energy accounting','Section coordinate x (m)','Energy (J)');
            hold(axEnergy,'on');
            plot(axEnergy,energy.x_m,energy.initial_kinetic_J+energy.gravity_work_J+energy.push_work_J, ...
                ':','Color',[0.45 0.30 0.70],'LineWidth',1.6);
            plot(axEnergy,energy.x_m,energy.rolling_loss_J+energy.viscous_loss_J, ...
                '--','Color',[0.90 0.35 0.08],'LineWidth',1.6);
            plot(axEnergy,energy.x_m,energy.viscous_loss_J,'-.', ...
                'Color',[0.65 0.18 0.48],'LineWidth',1.6);
            hold(axEnergy,'off');
            legend(axEnergy,{'Kinetic energy','Initial kinetic + gravity + push', ...
                'Coulomb + viscous losses','Viscous loss only'},'Location','best');
            plotResult(axAgreement,result.x_m,result.dynamic_minus_energy_mps, ...
                'Dynamic speed minus energy speed','Section coordinate x (m)','Speed difference (m/s)');
            hold(axAgreement,'on');
            plot(axAgreement,result.x_m,result.lagrangian_minus_dynamic_mps,'-.', ...
                'Color',[0.58 0.20 0.75],'LineWidth',1.8);
            hold(axAgreement,'off');
            title(axAgreement,'Speed differences (Lagrangian uses smoother terrain)');
            legend(axAgreement,{'Dynamic - energy','Lagrangian - dynamic'},'Location','best');
            plotResult(axLagrangianSpeed,lagrangian.time_s,lagrangian.speed_mps, ...
                'Lagrangian speed (physical, without display offset)','Time (s)','Speed (m/s)');
            plotResult(axLagrangianResidual,lagrangian.time_s,lagrangian.balance_residual_J, ...
                'K + potential change + losses - push work - K0','Time (s)','Energy residual (J)');
            lagDifference = result.lagrangian_minus_dynamic_mps;
            lagDifference = lagDifference(isfinite(lagDifference));
            if isempty(lagDifference), maxLagDifference = NaN; else, maxLagDifference = max(abs(lagDifference)); end
            lagrangianSummary.Text = sprintf(['%s Final: %.3f m/s at %.2f s.\n' ...
                'Max |Lagrangian - dynamic|: %.3g m/s | Max energy residual: %.3g J'], ...
                lagrangianMessage,lagrangian.speed_mps(end),lagrangian.time_s(end), ...
                maxLagDifference,max(abs(lagrangian.balance_residual_J)));
            difference = result.dynamic_minus_energy_mps;
            difference = difference(isfinite(difference));
            if isempty(difference), maxDifference = NaN; else, maxDifference = max(abs(difference)); end
            energySummary.Text = sprintf(['%s\n' ...
                'Max speed difference on shared x: %.3g m/s | Max energy residual: %.3g J'], ...
                energyMessage,maxDifference,max(abs(energy.balance_residual_J)));
            xHill = unique([linspace(p.xStart,p.xEnd,600),obstacleBreaks(p)]);
            elevationHill = hillHeight(xHill,p)+p.sourceElevationAtStart;
            plotResult(axHill,xHill,elevationHill, ...
                'Merged hill profile','Course coordinate x (m)','Elevation (m)');
            hold(axHill,'on');
            plot(axHill,result.x_m,result.original_elevation_m,'Color',[0.95 0.45 0.12],'LineWidth',2);
            plot(axHill,result.x_m(end),result.original_elevation_m(end),'o', ...
                'MarkerFaceColor',[0.95 0.45 0.12],'MarkerEdgeColor','w','MarkerSize',8);
            hold(axHill,'off');
            hold(axHill,'on');
            plot(axHill,xHill,hillHeight(xHill,pLag)+pLag.sourceElevationAtStart, ...
                '--','Color',[0.58 0.20 0.75],'LineWidth',1.6);
            hold(axHill,'off');
            legend(axHill,{'Dynamic / energy terrain','Dynamic path','Final position', ...
                'Lagrangian terrain'},'Location','best');
            plotResult(axGrade,xHill,-atan(hillSlope(xHill,p))*180/pi, ...
                'Local slope (positive downhill)','Section coordinate x (m)','Downhill angle (deg)');
            hold(axGrade,'on');
            plot(axGrade,xHill,-atan(hillSlope(xHill,pLag))*180/pi, ...
                '--','Color',[0.58 0.20 0.75],'LineWidth',1.6);
            hold(axGrade,'off');
            legend(axGrade,{'Dynamic / energy angle','Lagrangian angle (+30% smoothing)'},'Location','best');
            xlim(axEnergy,[p.xStart p.xEnd]);
            xlim(axAgreement,[p.xStart p.xEnd]);
            xlim(axHill,[p.xStart p.xEnd]);
            xlim(axGrade,[p.xStart p.xEnd]);
            markObstacles(axHill,p,true);
            metricSpeed.Text = sprintf('Final speed\n%.2f m/s | %.2f mph', ...
                result.displayed_speed_mps(end),result.displayed_speed_mps(end)*2.236936292);
            metricTime.Text = sprintf('Elapsed time\n%.2f s | %.2f m', ...
                result.time_s(end),result.distance_m(end));
            metricMass.Text = sprintf('Effective mass\n%.2f kg (+%.2f)', ...
                p.mEffective,p.mEffective-p.mass);
            updateCourseOverview(result,energy,lagrangian,p);
            prepareAnimation();
            statusLabel.Text = ['Merged hill: simulation ends at 116 m | ' message];
            exportButton.Enable = 'on';
            energyExportButton.Enable = 'on';
            lagrangianExportButton.Enable = 'on';
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
        ax.FontSize = 10; box(ax,'off');
    end

    function drawComparison(ax,result,energy,lagrangian,p,large)

        [dynamicX,dynamicValues] = trimComparisonSamples(result.comparison_x_m, ...
            [result.displayed_speed_mps result.dynamic_speed_lower_mps result.dynamic_speed_upper_mps]);
        [energyX,energySpeed] = trimComparisonSamples(energy.comparison_x_m,energy.displayed_speed_mps);
        [lagrangianX,lagrangianSpeed] = trimComparisonSamples( ...
            lagrangian.comparison_x_m,lagrangian.displayed_speed_mps);
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
        hLagrangian = plot(ax,lagrangianX,lagrangianSpeed,':', ...
            'Color',[0.58 0.20 0.75],'LineWidth',2.8,'Marker','none');
        hGPS = plot(ax,gpsX,gpsV,'-', ...
            'Color',gpsColor,'LineWidth',2.2,'Marker','none');
        hold(ax,'off');
        grid(ax,'on');

        xlim(ax,[min(0,gpsX(1)), ...
            max(p.xEnd+p.comparisonXOffset_m,gpsX(end))]);
        visibleValues = [dynamicValues(:,1);energySpeed;lagrangianSpeed;gpsV];
        if large
            visibleValues = [visibleValues;dynamicValues(:,2); ...
                dynamicValues(:,3);gpsLower;gpsUpper];
        end
        visibleValues = visibleValues(isfinite(visibleValues));
        top = max(visibleValues); bottom = min([0;visibleValues]);
        ylim(ax,[bottom-0.02*max(1,top-bottom),top+0.08*max(1,top-bottom)]);
        title(ax,'Simulated and recorded kart speed');
        xlabel(ax,'Comparison coordinate (m)'); ylabel(ax,'Speed (m/s)');
        labels = {'Dynamic','Energy balance','Real kart speed','Lagrangian'};
        if large
            ax.FontSize = 13;
            ax.UserData = [hDynamic hEnergy hGPS hBand hGPSBand hLagrangian];
            updateComparisonVisibility();
        else
            ax.FontSize = 10;
            legend(ax,[hDynamic hEnergy hGPS hLagrangian],labels,'Location','best','FontSize',8);
        end
    end

    function exportComparison(~,~)
        if ~strcmp(comparisonExportButton.Enable,'on'), return; end
        filters = {'*.csv','Selected curve and band data (*.csv)'; ...
            '*.png','Selected plot image (*.png)'; ...
            '*.pdf','Selected plot PDF (*.pdf)'};
        [name,folder,formatIndex] = uiputfile(filters,'Export speed comparison', ...
            'KartSpeedComparison_v48.csv');
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


                        rows = rows(isfinite(rows.x_m),:);
                        exported = [exported;rows];
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
        [name,folder] = uiputfile('*.csv','Export simulation results','KartVelocityResults_v48.csv');
        if isequal(name,0), return; end
        try

            exported = lastResults;
            keys = fieldnames(lastParameters);
            for j = 1:numel(keys)
                value = lastParameters.(keys{j});
                if isstruct(value)
                    exported.(['input_' keys{j} '_json']) = repmat({jsonencode(value)},height(exported),1);
                else
                    exported.(['input_' keys{j}]) = repmat(value,height(exported),1);
                end
            end
            writetable(exported,fullfile(folder,name));
            statusLabel.Text = ['Exported results and run parameters to ' fullfile(folder,name)];
        catch err
            uialert(fig,err.message,'Export error');
        end
    end

    function exportLagrangian(~,~)
        if isempty(lastLagrangian), return; end
        [name,folder] = uiputfile('*.csv','Export full Lagrangian results','KartLagrangianResults.csv');
        if isequal(name,0), return; end
        try
            exported = lastLagrangian;
            keys = fieldnames(lastLagrangianParameters);
            for j = 1:numel(keys)
                value = lastLagrangianParameters.(keys{j});
                if isstruct(value)
                    exported.(['input_' keys{j} '_json']) = repmat({jsonencode(value)},height(exported),1);
                else
                    exported.(['input_' keys{j}]) = repmat(value,height(exported),1);
                end
            end
            writetable(exported,fullfile(folder,name));
            statusLabel.Text = ['Exported Lagrangian results to ' fullfile(folder,name)];
        catch err
            uialert(fig,err.message,'Export error');
        end
    end

    function exportEnergy(~,~)
        if isempty(lastEnergy), return; end
        [name,folder] = uiputfile('*.csv','Export full energy results','KartEnergyResults_v48.csv');
        if isequal(name,0), return; end
        try
            exported = lastEnergy;
            keys = fieldnames(lastParameters);
            for j = 1:numel(keys)
                value = lastParameters.(keys{j});
                if isstruct(value)
                    exported.(['input_' keys{j} '_json']) = repmat({jsonencode(value)},height(exported),1);
                else
                    exported.(['input_' keys{j}]) = repmat(value,height(exported),1);
                end
            end
            writetable(exported,fullfile(folder,name));
            statusLabel.Text = ['Exported energy results to ' fullfile(folder,name)];
        catch err
            uialert(fig,err.message,'Export error');
        end
    end
end


function [results,p,status] = simulateKart(p)

p.sourceXStart = p.sourceXOffset+p.xStart;
p.sourceXEnd = p.sourceXOffset+p.xEnd;
[heightPP,~] = mergedHillInterpolants(p);
p.sourceElevationAtStart = ppval(heightPP,p.sourceXStart)+rampHeight(p.xStart,p);
p.mEffective = p.mass + p.nWheels*p.Iwheel/p.radius^2;
p.dragFactor = 0.5*p.rho*p.CdA;
p.linear1PerWheel_Nms = p.b1*p.radius^2/p.nWheels;
p.linear2PerWheel_Nms = p.b2*p.radius^2/p.nWheels;
p.linear3PerWheel_Nms = p.b3*p.radius^2/p.nWheels;
p.pushEndX_m = findPushEnd(p);


[initialAcceleration,~,~,~,~] = localForces(p.xStart,p.v0,p);
if p.v0 == 0 && initialAcceleration <= 0
    time_s = [0;p.tMax];
    state = [p.xStart 0 0;p.xStart 0 0];
    status = 'Remained at rest: local downhill gravity plus push does not exceed assumed rolling resistance. No rollback is modeled.';
else
    options = odeset('RelTol',1e-8,'AbsTol',[1e-9 1e-10 1e-9],'MaxStep',0.05);
    breaks = obstacleBreaks(p);
    time_s = 0;
    state = [p.xStart p.v0 0];
    eventIndex = [];
    reachedEnd = false;
    for segment = 1:numel(breaks)-1
        segmentParameters = p;
        segmentParameters.segmentMid = (breaks(segment)+breaks(segment+1))/2;
        segmentOptions = odeset(options,'Events', ...
            @(t,y) stopEvents(t,y,p,breaks(segment+1)));
        [tt,yy,~,~,eventIndex] = ode45(@(t,y) kartODE(t,y,segmentParameters), ...
            [time_s(end) p.tMax],state(end,:).',segmentOptions);
        time_s = [time_s;tt(2:end)];
        state = [state;yy(2:end,:)];
        if any(eventIndex == 1), state(end,1) = breaks(segment+1); end
        reachedEnd = any(eventIndex == 1) && segment == numel(breaks)-1;
        if any(eventIndex == 2) || ~any(eventIndex == 1) || time_s(end) >= p.tMax
            break
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
original_x_m = x_m+p.sourceXOffset;
horizontal_from_start_m = x_m-p.xStart;
original_elevation_m = elevation_m+p.sourceElevationAtStart;
speed_mps = max(state(:,2),0);
speed_mph = speed_mps*2.236936292;
[acceleration_mps2,downhill_angle_deg,gravity_force_N,rolling_force_N,drag_force_N,viscous_force_N,push_force_N] = ...
    localForces(x_m,speed_mps,p);
acceleration_mps2(speed_mps <= 0 & acceleration_mps2 < 0) = 0;
results = table(time_s,x_m,original_x_m,horizontal_from_start_m,distance_m,elevation_m,original_elevation_m,speed_mps,speed_mph, ...
    acceleration_mps2,downhill_angle_deg,gravity_force_N,rolling_force_N,drag_force_N,viscous_force_N,push_force_N);
results.push_work_J = p.pushForce_N*min(distance_m,p.pushDistance_m);
results.mu_coefficient = surfaceMu(x_m,p);
[~,results.viscous_b_Ns_per_m,results.friction_section] = sectionCoefficients(x_m,p);
end

function dydt = kartODE(~,y,p)
v = max(y(2),0);
slope = hillSlope(y(1),p);
a = localForces(y(1),v,p);

dydt = [v/sqrt(1+slope^2);a;v];
end

function [a,angleDeg,Fgravity,Frolling,Fdrag,Fviscous,Fpush] = localForces(x,v,p)
slope = hillSlope(x,p);
q = sqrt(1+slope.^2);
angleDeg = -atan(slope)*180/pi;
Fgravity = -p.mass*9.81*slope./q;
[Frolling,Fviscous] = frictionForces(v,q,p,x);
Fdrag = p.dragFactor*v.^2;

Fpush = startPushForce(x,p);
a = (Fgravity+Fpush-Frolling-Fdrag-Fviscous)/p.mEffective;
end

function [Frolling,Fviscous] = frictionForces(v,q,p,x)
speed = max(v,0);
Frolling = surfaceMu(x,p).*p.mass*9.81./q;
[~,b] = sectionCoefficients(x,p);
Fviscous = b.*speed;
end

function [value,isterminal,direction] = stopEvents(t,y,p,xTarget)
stopValue = y(2);


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


function [results,status] = simulateLagrangian(p)


q0 = sqrt(1+hillSlope(p.xStart,p)^2);
y0 = [p.xStart;p.v0/q0;0;0;0;0;0];
initialDerivative = lagrangianODE(0,y0,p);
if p.v0 == 0 && initialDerivative(2) <= 0
    time_s = [0;p.tMax];
    state = [y0.';y0.'];
    status = 'Lagrangian: remained at rest.';
else
    options = odeset('RelTol',1e-10, ...
        'AbsTol',[1e-10 1e-12 1e-10 1e-9 1e-9 1e-9 1e-9],'MaxStep',0.02);


    [~,terrainSlopePP] = mergedHillInterpolants(p);
    terrainKnots = terrainSlopePP.breaks-p.sourceXOffset;
    terrainKnots = terrainKnots(terrainKnots>p.xStart & terrainKnots<p.xEnd);
    breaks = unique([obstacleBreaks(p) terrainKnots]);
    time_s = 0;
    state = y0.';
    eventIndex = [];
    reachedEnd = false;
    for segment = 1:numel(breaks)-1
        segmentParameters = p;
        segmentParameters.segmentMid = (breaks(segment)+breaks(segment+1))/2;
        segmentOptions = odeset(options,'Events', ...
            @(t,y) stopEvents(t,y,p,breaks(segment+1)));
        [tt,yy,~,~,eventIndex] = ode45(@(t,y) lagrangianODE(t,y,segmentParameters), ...
            [time_s(end) p.tMax],state(end,:).',segmentOptions);
        time_s = [time_s;tt(2:end)];
        state = [state;yy(2:end,:)];
        if any(eventIndex == 1)
            state(end,1) = breaks(segment+1);


            qLeft = sqrt(1+hillSlope(state(end,1),segmentParameters)^2);
            qRight = sqrt(1+hillSlope(state(end,1),p)^2);
            state(end,2) = state(end,2)*qLeft/qRight;
        end
        reachedEnd = any(eventIndex == 1) && segment == numel(breaks)-1;
        if any(eventIndex == 2) || ~any(eventIndex == 1) || time_s(end) >= p.tMax
            break
        end
    end
    if reachedEnd
        status = 'Lagrangian: reached end x.';
    elseif any(eventIndex == 2)
        state(end,2) = 0;
        status = 'Lagrangian: stopped before end x (no rollback).';
    else
        status = 'Lagrangian: reached time limit.';
    end
end
x_m = state(:,1);
horizontal_velocity_mps = max(state(:,2),0);
distance_m = state(:,3);
slope = hillSlope(x_m,p);
curvature = hillSecondDerivative(x_m,p);
q = sqrt(1+slope.^2);
speed_mps = q.*horizontal_velocity_mps;
[rolling_force_N,viscous_force_N] = frictionForces(speed_mps,q,p,x_m);
drag_force_N = p.dragFactor*speed_mps.^2;
push_force_N = startPushForce(x_m,p);
generalized_force_N = (push_force_N-rolling_force_N-viscous_force_N-drag_force_N).*q;
horizontal_acceleration_mps2 = (generalized_force_N-p.mass*9.81*slope ...
    -p.mEffective*slope.*curvature.*horizontal_velocity_mps.^2)./(p.mEffective*q.^2);
acceleration_mps2 = q.*horizontal_acceleration_mps2 ...
    +slope.*curvature./q.*horizontal_velocity_mps.^2;
stationary = speed_mps <= 0 & acceleration_mps2 < 0;
horizontal_acceleration_mps2(stationary) = 0;
acceleration_mps2(stationary) = 0;
elevation_m = hillHeight(x_m,p);
original_x_m = x_m+p.sourceXOffset;
kinetic_J = 0.5*p.mEffective*speed_mps.^2;
initial_kinetic_J = repmat(0.5*p.mEffective*p.v0^2,numel(x_m),1);
potential_change_J = p.mass*9.81*(elevation_m-hillHeight(p.xStart,p));
rolling_loss_J = state(:,4);
viscous_loss_J = state(:,5);
drag_loss_J = state(:,6);
push_work_J = state(:,7);
balance_residual_J = kinetic_J+potential_change_J+rolling_loss_J ...
    +viscous_loss_J+drag_loss_J-push_work_J-initial_kinetic_J;
results = table(time_s,x_m,original_x_m,distance_m,elevation_m, ...
    horizontal_velocity_mps,horizontal_acceleration_mps2,speed_mps,acceleration_mps2, ...
    rolling_force_N,viscous_force_N,drag_force_N,generalized_force_N, ...
    kinetic_J,initial_kinetic_J,potential_change_J,rolling_loss_J,viscous_loss_J, ...
    drag_loss_J,push_force_N,push_work_J,balance_residual_J);
results.mu_coefficient = surfaceMu(x_m,p);
[~,results.viscous_b_Ns_per_m,results.friction_section] = sectionCoefficients(x_m,p);
end

function dydt = lagrangianODE(~,y,p)
slope = hillSlope(y(1),p);
curvature = hillSecondDerivative(y(1),p);
q = sqrt(1+slope^2);
u = max(y(2),0);
v = q*u;
[Frolling,Fviscous] = frictionForces(v,q,p,y(1));
Fdrag = p.dragFactor*v^2;
Fpush = startPushForce(y(1),p);
Qx = (Fpush-Frolling-Fviscous-Fdrag)*q;
udot = (Qx-p.mass*9.81*slope-p.mEffective*slope*curvature*u^2) ...
    /(p.mEffective*q^2);
dydt = [u;udot;v;Frolling*v;Fviscous*v;Fdrag*v;Fpush*v];
end

function curvature = hillSecondDerivative(x,p)


[heightPP,~,curvaturePP] = mergedHillInterpolants(p);
xSource = min(max(x+p.sourceXOffset,heightPP.breaks(1)),heightPP.breaks(end));
curvature = ppval(curvaturePP,xSource);
end


function [results,status] = simulateEnergy(p)


K0 = 0.5*p.mEffective*p.v0^2;
initialSlope = hillSlope(p.xStart,p);
initialGradient = -p.mass*9.81*(initialSlope+surfaceMu(p.xStart,p)) ...
    + startPushForce(p.xStart,p)*sqrt(1+initialSlope^2);
if K0 == 0 && initialGradient <= 0
    x_m = p.xStart;
    state = [0 0 0 0];
    status = 'Energy method: remained at rest.';
else
    options = odeset('RelTol',1e-9,'AbsTol',[1e-7 1e-7 1e-7 1e-7], ...
        'MaxStep',min(0.05,(p.xEnd-p.xStart)/1000), ...
        'Events',@(x,z) energyStop(x,z,p,K0));
    breaks = obstacleBreaks(p);
    x_m = p.xStart;
    state = [K0 0 0 0];
    eventIndex = [];
    for segment = 1:numel(breaks)-1
        segmentParameters = p;
        segmentParameters.segmentMid = (breaks(segment)+breaks(segment+1))/2;
        [xx,zz,~,~,eventIndex] = ode45(@(x,z) energyODE(x,z,segmentParameters), ...
            [breaks(segment) breaks(segment+1)],state(end,:).',options);
        x_m = [x_m;xx(2:end)];
        state = [state;zz(2:end,:)];
        if ~isempty(eventIndex), break; end
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
original_x_m = x_m+p.sourceXOffset;
horizontal_from_start_m = x_m-p.xStart;
original_elevation_m = elevation_m+p.sourceElevationAtStart;
initial_kinetic_J = repmat(K0,numel(x_m),1);
gravity_work_J = p.mass*9.81*(hillHeight(p.xStart,p)-elevation_m);
rolling_loss_J = coulombWork(x_m,p);
drag_loss_J = state(:,2);
viscous_loss_J = state(:,3);
push_work_J = state(:,4);
push_force_N = startPushForce(x_m,p);

balance_residual_J = kinetic_J-initial_kinetic_J-gravity_work_J-push_work_J+rolling_loss_J+drag_loss_J+viscous_loss_J;
results = table(x_m,original_x_m,horizontal_from_start_m,elevation_m,original_elevation_m,speed_mps, ...
    kinetic_J,initial_kinetic_J,gravity_work_J,rolling_loss_J,drag_loss_J,viscous_loss_J,push_work_J,push_force_N,balance_residual_J);
results.mu_coefficient = surfaceMu(x_m,p);
[~,results.viscous_b_Ns_per_m,results.friction_section] = sectionCoefficients(x_m,p);
end

function dzdx = energyODE(x,z,p)
slope = hillSlope(x,p);
q = sqrt(1+slope^2);
speedSquared = 2*max(z(1),0)/p.mEffective;
dragWorkGradient = p.dragFactor*speedSquared*q;
[Frolling,Fviscous] = frictionForces(sqrt(speedSquared),q,p,x);
viscousWorkGradient = Fviscous*q;
pushWorkGradient = startPushForce(x,p)*q;
dzdx = [-p.mass*9.81*slope-Frolling*q-dragWorkGradient-viscousWorkGradient+pushWorkGradient; ...
    dragWorkGradient;viscousWorkGradient;pushWorkGradient];
end

function [value,isterminal,direction] = energyStop(x,z,p,K0)
value = z(1);
if x == p.xStart && K0 == 0, value = 1; end
isterminal = 1;
direction = -1;
end

function sampled = sampleAtX(x,values,query)


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


function h = hillHeight(x,p)
[heightPP,~] = mergedHillInterpolants(p);

xSource = min(max(x+p.sourceXOffset,heightPP.breaks(1)),heightPP.breaks(end));
h = ppval(heightPP,xSource)+rampHeight(x,p)-p.sourceElevationAtStart;
end

function slope = hillSlope(x,p)
[heightPP,slopePP] = mergedHillInterpolants(p);
xSource = min(max(x+p.sourceXOffset,heightPP.breaks(1)),heightPP.breaks(end));
slope = ppval(slopePP,xSource)+rampSlope(x,p);
end

function [heightPP,slopePP,curvaturePP] = mergedHillInterpolants(p)
persistent heightCache slopeCache curvatureCache smoothingKey
key = [p.angleSmoothingFraction p.angleSmoothingSigma_m];
if isempty(heightCache) || ~isequal(smoothingKey,key)


    data = [
        0 8;
        0.24995715627450296 7.9203675976386885;
        0.49991431254900592 7.8426462867709308;
        0.74987146882350886 7.7666505411839815;
        0.99982862509801185 7.6922764017868825;
        1.2497857813725148 7.6194793460544457;
        1.4997429376470177 7.5482536521215433;
        1.7497000939215208 7.4786141918770781;
        1.9996572501960237 7.4105813885466363;
        2.2496144064705268 7.3441698058720997;
        2.4995715627450297 7.2793805350189711;
        2.7495287190195326 7.2161972494570978;
        2.9994858752940354 7.1545855402636356;
        3.2494430315685383 7.0944949487288742;
        3.4994001878430416 7.0358629928270897;
        3.7493573441175445 6.9786204244142702;
        3.9993145003920474 6.9226970328063704;
        4.2492716566665507 6.8680273416268962;
        4.4992288129410536 6.8145556359959283;
        4.7491859692155565 6.7622399667827233;
        4.9991431254900593 6.7110548664160712;
        5.2491002817645622 6.6609926501229504;
        5.4990574380390651 6.6120632896811742;
        5.749014594313568 6.5642929419225426;
        5.9989717505880709 6.5177212912512363;
        6.2489289068625737 6.4723979256400392;
        6.4988860631370766 6.4283780103231445;
        6.7488432194115804 6.3857175533755806;
        6.9988003756860833 6.3444685722072052;
        7.2487575319605861 6.3046744686066623;
        7.498714688235089 6.2663659011475001;
        7.7486718445095919 6.2295574070188602;
        7.9986290007840948 6.1942449716074544;
        8.2485861570585985 6.1604046762632887;
        8.4985433133331014 6.1279949125701316;
        8.7485004696076043 6.0969510135065876;
        8.9984576258821072 6.0671903415113206;
        9.2484147821566101 6.0386179794501382;
        9.4983719384311129 6.0111275554234913;
        9.7483290947056158 5.9846049888341568;
        9.9982862509801187 5.9589321221489966;
        10.248243407254622 5.9339907214258538;
        10.498200563529124 5.9096657034861764;
        10.748157719803627 5.8858479863011954;
        10.99811487607813 5.8624367323840794;
        11.248072032352633 5.8393409286040967;
        11.498029188627136 5.8164802961239479;
        11.747986344901639 5.793785565691163;
        11.997943501176142 5.7711981853154546;
        12.247900657450645 5.7486695486795227;
        12.497857813725147 5.7261598440489578;
        12.74781496999965 5.7036366485477439;
        12.997772126274153 5.6810733908582174;
        13.247729282548656 5.6584475460730568;
        13.497686438823161 5.6357390222017099;
        13.747643595097664 5.6129285708264751;
        13.997600751372167 5.5899963500012735;
        14.247557907646669 5.566920702213964;
        14.497515063921172 5.5436773221395246;
        14.747472220195675 5.5202396483643614;
        14.997429376470178 5.4965770704035934;
        15.247386532744681 5.4726559355215407;
        15.497343689019184 5.4484403259747465;
        15.747300845293687 5.4238933387243407;
        15.99725800156819 5.3989786254157828;
        16.247215157842692 5.3736616381869711;
        16.497172314117197 5.3479141797846896;
        16.747129470391698 5.321716440217096;
        16.997086626666203 5.295060013743794;
        17.247043782940704 5.2679510852602069;
        17.497000939215209 5.2404135139217134;
        17.74695809548971 5.2124916321292636;
        17.996915251764214 5.1842525621673863;
        18.246872408038715 5.1557878428454869;
        18.49682956431322 5.1272141557402495;
        18.746786720587721 5.0986729470607708;
        18.996743876862226 5.0703287595853741;
        19.246701033136727 5.0423661230725294;
        19.496658189411232 5.014984904729018;
        19.746615345685733 4.9883940967117866;
        19.996572501960237 4.9628041261677831;
        20.246529658234738 4.9384178509914385;
        20.496486814509243 4.9154206762382939;
        20.746443970783744 4.893970661307419;
        20.996401127058249 4.8741873605280945;
        21.24635828333275 4.8561427693941903;
        21.496315439607255 4.8398534752612745;
        21.746272595881759 4.825275148177302;
        21.99622975215626 4.8123000344883513;
        22.246186908430765 4.8007579527408684;
        22.496144064705266 4.790421058284644;
        22.746101220979771 4.7810123579877057;
        22.996058377254272 4.7722176467740693;
        23.246015533528777 4.7637002363245413;
        23.495972689803278 4.7551175633355314;
        23.745929846077782 4.7461385547840713;
        23.995887002352283 4.7364609210194555;
        24.245844158626788 4.7258267121276942;
        24.495801314901289 4.7140354464726348;
        24.745758471175794 4.7009539180138376;
        24.995715627450295 4.6865221408180862;
        25.2456727837248 4.6707551929266469;
        25.495629939999301 4.6537410234007028;
        25.745587096273805 4.6356345565475694;
        25.995544252548306 4.6166486439529999;
        26.245501408822811 4.5970427398202283;
        26.495458565097312 4.5771094239598824;
        26.745415721371817 4.5571605129585926;
        26.995372877646322 4.5375129673615051;
        27.245330033920823 4.5184753582928021;
        27.495287190195327 4.5003355285243547;
        27.745244346469828 4.4833499339077125;
        27.995201502744333 4.4677350359533294;
        28.245158659018834 4.4536610005421506;
        28.495115815293339 4.4412478431388305;
        28.74507297156784 4.430564048007632;
        28.995030127842345 4.4216275652718906;
        29.244987284116846 4.414408996665915;
        29.49494444039135 4.4088368432090519;
        29.744901596665851 4.4048040849839332;
        29.994858752940356 4.4021759616746179;
        30.244815909214857 4.4007983790191947;
        30.494773065489362 4.400506451724457;
        30.744730221763863 4.401132641490145;
        30.994687378038368 4.4025136712325033;
        31.244644534312869 4.4044984802588099;
        31.494601690587373 4.4069521157805127;
        31.744558846861874 4.4097592796085339;
        31.994516003136379 4.4128263310157765;
        32.24447315941088 4.4160819648934098;
        32.494430315685385 4.4194767395929757;
        32.74438747195989 4.4229797574413512;
        32.994344628234394 4.4265770045000625;
        33.244301784508892 4.4302680345503562;
        33.494258940783396 4.4340624162183158;
        33.744216097057901 4.437976154756119;
        33.994173253332406 4.4420284466076305;
        34.244130409606903 4.4462388041583569;
        34.494087565881408 4.450624696221773;
        34.744044722155913 4.4551997953468838;
        34.994001878430417 4.4599728276890369;
        35.243959034704915 4.4649469888306452;
        35.493916190979419 4.4701198533973452;
        35.743873347253924 4.4754836822137785;
        35.993830503528429 4.4810260184175661;
        36.243787659802926 4.4867304626732638;
        36.493744816077431 4.4925775258078327;
        36.743701972351936 4.4985454236571609;
        36.99365912862644 4.5046100798680708;
        37.243616284900938 4.510747309239485;
        37.493573441175442 4.5169323628027058;
        37.743530597449947 4.5231402319177239;
        37.993487753724452 4.5293457790027079;
        38.243444909998956 4.5355237131632169;
        38.493402066273454 4.5416484886209956;
        38.743359222547959 4.547694537957681;
        38.993316378822463 4.5536351697649931;
        39.243273535096968 4.5594425651357291;
        39.493230691371465 4.5650873987322882;
        39.74318784764597 4.5705384641113023;
        39.993145003920475 4.5757623373719305;
        40.243102160194979 4.5807228637091004;
        40.493059316469477 4.5853810008366596;
        40.743016472743982 4.5896946456713001;
        40.992973629018486 4.5936186327370994;
        41.242930785292991 4.5971049744160695;
        41.492887941567488 4.6001034218949854;
        41.742845097841993 4.6025624265816436;
        41.992802254116498 4.6044305713645484;
        42.242759410391002 4.6056585166757591;
        42.4927165666655 4.6062014270621372;
        42.742673722940005 4.6060219314050288;
        42.992630879214509 4.6050934908360208;
        43.242588035489014 4.6034038006493381;
        43.492545191763519 4.600958200521867;
        43.742502348038016 4.5977827245732703;
        43.992459504312521 4.5939264734608738;
        44.242416660587025 4.5894629973282974;
        44.49237381686153 4.5844904215225011;
        44.742330973136028 4.5791301265432107;
        44.992288129410532 4.5735239044859952;
        45.242245285685037 4.5678296460793746;
        45.492202441959542 4.56221562281298;
        45.742159598234039 4.5568539935576258;
        45.992116754508544 4.5519138084077397;
        46.242073910783049 4.5475537840396676;
        46.492031067057553 4.5439156058059442;
        46.741988223332051 4.5411181484277892;
        46.991945379606555 4.5392529993652673;
        47.24190253588106 4.5383815490690775;
        47.491859692155565 4.5385337698016865;
        47.741816848430062 4.5397086565446019;
        47.991774004704567 4.5418761701334533;
        48.241731160979072 4.5449803978551229;
        48.491688317253576 4.5489435873798678;
        48.741645473528081 4.553670688274817;
        48.991602629802578 4.5590540069186236;
        49.241559786077083 4.5649776605555168;
        49.491516942351588 4.5713215594722607;
        49.741474098626092 4.5779647518354309;
        49.99143125490059 4.5847880989295255;
        50.241388411175095 4.5916761627786755;
        50.491345567449599 4.5985185420891721;
        50.741302723724104 4.605210725521431;
        50.991259879998601 4.6116545927020827;
        51.241217036273106 4.6177587445951707;
        51.491174192547611 4.6234386469721596;
        51.741131348822115 4.6286171874267232;
        51.991088505096613 4.6332251368237456;
        52.241045661371118 4.6372017706017106;
        52.491002817645622 4.640495623424278;
        52.740959973920127 4.6430652779046415;
        52.990917130194624 4.6448800851531669;
        53.240874286469129 4.6459209979358311;
        53.490831442743634 4.6461800939176472;
        53.740788599018138 4.645660894729323;
        53.990745755292643 4.6443778253680241;
        54.240702911567141 4.6423550991949112;
        54.490660067841645 4.6396251224196403;
        54.74061722411615 4.6362265964267158;
        54.990574380390655 4.6322020501545902;
        55.240531536665152 4.6275962601129619;
        55.490488692939657 4.6224533322302843;
        55.740445849214161 4.6168144848490051;
        55.990403005488666 4.6107161073709371;
        56.240360161763164 4.604188227057187;
        56.490317318037668 4.5972534897628634;
        56.740274474312173 4.5899267267105222;
        56.990231630586678 4.5822151409315612;
        57.240188786861175 4.5741191059477568;
        57.49014594313568 4.5656335278780897;
        57.740103099410184 4.5567497328960567;
        57.990060255684689 4.5474576006195964;
        58.240017411959187 4.5377479732075736;
        58.489974568233691 4.5276150933566406;
        58.739931724508196 4.5170588983297231;
        58.989888880782701 4.5060869891947801;
        59.239846037057205 4.4947161112000344;
        59.489803193331703 4.4829730926661311;
        59.739760349606208 4.4708950948010378;
        59.989717505880712 4.4585291699802418;
        60.239674662155217 4.4459311562095412;
        60.489631818429714 4.4331639884198477;
        60.739588974704219 4.4202955534582129;
        60.989546130978724 4.4073962504571389;
        61.239503287253228 4.3945362992641925;
        61.489460443527726 4.3817833765306204;
        61.739417599802231 4.3692004504604061;
        61.989374756076735 4.356843735163423;
        62.23933191235124 4.3447612081590288;
        62.489289068625737 4.3329918151312672;
        62.739246224900242 4.3215650666447152;
        62.989203381174747 4.310501100577353;
        63.239160537449251 4.2998111224605466;
        63.489117693723749 4.2894981110031498;
        63.739074849998254 4.2795576653569229;
        63.989032006272758 4.2699788744083271;
        64.238989162547256 4.2607451051660181;
        64.48894631882176 4.2518346970058403;
        64.738903475096265 4.2432214473504439;
        64.98886063137077 4.2348746962060417;
        65.238817787645274 4.2267594894458345;
        65.488774943919779 4.2188365486377659;
        65.738732100194284 4.2110622262003554;
        65.988689256468788 4.2033885319394901;
        66.238646412743279 4.1957633153808125;
        66.488603569017783 4.1881307086876474;
        66.738560725292288 4.1804318892839669;
        66.988517881566793 4.1726062011431102;
        67.238475037841297 4.1645926404533578;
        67.488432194115802 4.1563316752635204;
        67.738389350390307 4.1477674280019396;
        67.988346506664811 4.1388496958868188;
        68.238303662939302 4.1295362596854588;
        68.488260819213806 4.1197952605083614;
        68.738217975488311 4.1096070477914308;
        68.988175131762816 4.0989656979053315;
        69.23813228803732 4.0878799962834567;
        69.488089444311825 4.0763737776395761;
        69.73804660058633 4.0644855556225483;
        69.988003756860834 4.0522674167946189;
        70.237960913135339 4.0397831537382851;
        70.48791806940983 4.02710586567175;
        70.737875225684334 4.0143150049948222;
        70.987832381958839 4.0014929245892352;
        71.237789538233343 3.9887212964002057;
        71.487746694507848 3.9760775129752379;
        71.737703850782353 3.9636312668709857;
        71.987661007056857 3.9514414927454959;
        72.237618163331362 3.9395538947075788;
        72.487575319605853 3.9279990835374048;
        72.737532475880357 3.9167913641355638;
        72.987489632154862 3.9059284822216824;
        73.237446788429367 3.8953921196349595;
        73.487403944703871 3.8851491535397864;
        73.737361100978376 3.8751535977543266;
        73.98731825725288 3.8653490802105446;
        74.237275413527385 3.8556718403814001;
        74.487232569801876 3.8460538989390689;
        74.73718972607638 3.8364263480908778;
        74.987146882350885 3.8267225928602988;
        75.23710403862539 3.8168813803043609;
        75.487061194899894 3.8068494849743311;
        75.737018351174399 3.7965839434619761;
        75.986975507448903 3.786053759257781;
        76.236932663723408 3.7752410293738228;
        76.486889819997913 3.7641414743060349;
        76.736846976272403 3.7527643811897025;
        76.986804132546908 3.7411318898915664;
        77.236761288821413 3.7292779592400662;
        77.486718445095917 3.7172469434810198;
        77.736675601370422 3.7050914925235849;
        77.986632757644927 3.6928704177332632;
        78.236589913919431 3.6806464210848389;
        78.486547070193936 3.6684836723161975;
        78.736504226468426 3.6564453916310984;
        78.986461382742931 3.6445915191888947;
        79.236418539017436 3.6329765465469785;
        79.48637569529194 3.6216475774687349;
        79.736332851566445 3.6106426760128199;
        79.98629000784095 3.5999895485438831;
        80.236247164115454 3.5897045933887211;
        80.486204320389959 3.5797923377177021;
        80.736161476664464 3.5702452288835476;
        80.986118632938954 3.5610438536875848;
        81.236075789213459 3.5521576522679386;
        81.486032945487963 3.5435456087826993;
        81.735990101762468 3.5351575398216588;
        81.985947258036973 3.5269354929770409;
        82.235904414311477 3.5188153473376906;
        82.485861570585982 3.5107285629786196;
        82.735818726860487 3.5026040192472485;
        82.985775883134977 3.4943698911771395;
        83.235733039409482 3.4859555117107286;
        83.485690195683986 3.4772931648321439;
        83.735647351958491 3.4683197517802147;
        83.985604508232996 3.4589782703926399;
        84.2355616645075 3.4492190659391744;
        84.485518820782005 3.4390007924114472;
        84.73547597705651 3.4282908442509687;
        84.985433133331 3.4170656740925578;
        85.235390289605505 3.4053105994619726;
        85.485347445880009 3.3930192811362323;
        85.735304602154514 3.3801929698140172;
        85.985261758429019 3.3668390715303667;
        86.235218914703523 3.3529698634613863;
        86.485176070978028 3.3386012256039335;
        86.735133227252533 3.3237512259738415;
        86.985090383527037 3.3084389075949496;
        87.235047539801528 3.2926832742997147;
        87.485004696076032 3.2765026934761496;
        87.734961852350537 3.259914562739278;
        87.984919008625042 3.2429353683389737;
        88.234876164899546 3.2255810821109203;
        88.484833321174051 3.2078678216325414;
        88.734790477448556 3.1898127794960343;
        88.98474763372306 3.1714353229240047;
        89.234704789997551 3.1527580073624084;
        89.484661946272055 3.1338074199541412;
        89.73461910254656 3.1146149786311126;
        89.984576258821065 3.0952174154823791;
        90.234533415095569 3.0756569895792825;
        90.484490571370074 3.0559811673617969;
        90.734447727644579 3.0362420581308323;
        90.984404883919083 3.0164954398560337;
        91.234362040193588 2.9967993146279053;
        91.484319196468078 2.9772125677539751;
        91.734276352742583 2.9577931458078663;
        91.984233509017088 2.9385962613350705;
        92.234190665291592 2.9196725043953751;
        92.484147821566097 2.9010664693444164;
        92.734104977840602 2.8828153444643321;
        92.984062134115106 2.8649479432872575;
        93.234019290389611 2.8474840164164359;
        93.483976446664101 2.8304343032481181;
        93.733933602938606 2.8138006874057524;
        93.983890759213111 2.7975768985299543;
        94.233847915487615 2.7817492014199945;
        94.48380507176212 2.7662977251357774;
        94.733762228036625 2.7511973546654351;
        94.983719384311129 2.7364186094282532;
        95.233676540585634 2.721928492205564;
        95.483633696860124 2.7076903978575531;
        95.733590853134629 2.6936641208332865;
        95.983548009409134 2.6798053130981554;
        96.233505165683638 2.6660646303820856;
        96.483462321958143 2.652386563479348;
        96.733419478232648 2.6387082512166233;
        96.983376634507152 2.6249583659351874;
        97.233333790781657 2.6110566191559279;
        97.483290947056162 2.5969128869337053;
        97.733248103330652 2.5824277240786491;
        97.983205259605157 2.5674932315995451;
        98.233162415879661 2.5519944268760244;
        98.483119572154166 2.5358124361388006;
        98.733076728428671 2.51882745266012;
        98.983033884703175 2.5009225619419064;
        99.23299104097768 2.4819878964519808;
        99.482948197252185 2.4619252032208738;
        99.732905353526675 2.440652177664024;
        99.98286250980118 2.4181065385714691;
        100.23281966607568 2.3942493121946749;
        100.48277682235019 2.3690678164340473;
        100.73273397862469 2.3425771471862951;
        100.9826911348992 2.3148205792028569;
        101.2326482911737 2.2858684336725261;
        101.48260544744821 2.2558160487685428;
        101.73256260372271 2.2247799482489796;
        101.9825197599972 2.1928928563708174;
        102.23247691627171 2.1602982094051741;
        102.48243407254621 2.1271429269660302;
        102.73239122882072 2.0935702184897131;
        102.98234838509522 2.059712389615973;
        103.23230554136973 2.0256841530025991;
        103.48226269764423 1.9915760350257572;
        103.73221985391874 1.9574493707625755;
        103.98217701019323 1.923332672829182;
        104.23213416646773 1.8892197109329132;
        104.48209132274224 1.8550692643785303;
        104.73204847901674 1.8208068990974935;
        104.98200563529124 1.7863285329271221;
        105.23196279156575 1.7515057653062345;
        105.48191994784025 1.7161922493718949;
        105.73187710411476 1.6802314937380727;
        105.98183426038925 1.6434650465273015;
        106.23179141666375 1.6057406190051624;
        106.48174857293826 1.5669207257579847;
        106.73170572921276 1.526890129188909;
        106.98166288548727 1.4855623128432935;
        107.23162004176177 1.4428855664252445;
        107.48157719803628 1.3988455580556636;
        107.73153435431078 1.3534677068488017;
        107.98149151058529 1.306817661830886;
        108.23144866685978 1.2589995556939149;
        108.48140582313428 1.2101535982487246;
        108.73136297940879 1.1604511690231585;
        108.98132013568329 1.1100892745894539;
        109.2312772919578 1.0592838031268532;
        109.4812344482323 1.0082626689590029;
        109.7311916045068 0.95725804709829421;
        109.98114876078131 0.90649935563186135;
        110.2311059170558 0.85620694841518918;
        110.4810630733303 0.80658594698089814;
        110.73102022960481 0.75782154726317863;
        110.98097738587931 0.71007595856305028;
        111.23093454215382 0.66348731406013184;
        111.48089169842832 0.61816808284947022;
        111.73084885470283 0.57420683562654018;
        111.98080601097733 0.53167007866188609;
        112.23076316725184 0.490604490474287;
        112.48072032352633 0.45104121299438876;
        112.73067747980083 0.41299900280434626;
        112.98063463607534 0.3764875590958604;
        113.23059179234984 0.34151050304980107;
        113.48054894862435 0.30806762208227589;
        113.73050610489885 0.27615725102002164;
        113.98046326117336 0.24577726751435883;
        114.23042041744786 0.21692601341776213;
        114.48037757372235 0.18960321438034961;
        114.73033472999686 0.16380987764278399;
        114.98029188627136 0.139549320321711;
        115.23024904254586 0.11682786772078256;
        115.48020619882037 0.095655714673171965;
        115.73016335509487 0.076047607408277804;
        115.98012051136938 0.058023560810383283;
        116.23007766764388 0.041609077774430681;
        116.48003482391837 0.026834855883481976;
        116.72999198019288 0.013735215345387727;
        116.97994913646738 0.0023461349604643623;
        117.22990629274189 0;
        117.47986344901639 0;
        117.7298206052909 0;
        117.9797777615654 0;
        118.22973491783991 0;
        118.47969207411441 0;
        118.7296492303889 0;
        118.97960638666341 0;
        119.22956354293791 0;
        119.47952069921242 0;
        119.72947785548692 0;
        119.97943501176142 0;
    ];


    rawHeight = pchip(data(:,1),data(:,2));
    [breaks,coefs,~,order,dimension] = unmkpp(rawHeight);
    rawSlope = mkpp(breaks,bsxfun(@times,coefs(:,1:order-1),order-1:-1:1),dimension);
    slopeSamples = ppval(rawSlope,data(:,1));
    slopeSamples = slopeSamples(:);
    sigma_m = 0.75;
    spacing = median(diff(data(:,1)));
    radius = ceil(3*sigma_m/spacing);
    offsets = (-radius:radius)'*spacing;
    weights = exp(-0.5*(offsets/sigma_m).^2);
    weights = weights/sum(weights);
    padded = [repmat(slopeSamples(1),radius,1);slopeSamples; ...
        repmat(slopeSamples(end),radius,1)];
    smoothSlope = conv(padded,weights,'valid');


    for pass = 1:(1+double(p.angleSmoothingFraction>0))
        if pass == 2
            angles = -atan(ppval(slopeCache,data(:,1)));
            angles = angles(:);
            radius = ceil(3*p.angleSmoothingSigma_m/spacing);
            offsets = (-radius:radius)'*spacing;
            weights = exp(-0.5*(offsets/p.angleSmoothingSigma_m).^2);
            weights = weights/sum(weights);
            paddedAngles = [repmat(angles(1),radius,1);angles;repmat(angles(end),radius,1)];
            filteredAngles = conv(paddedAngles,weights,'valid');
            angles = angles+p.angleSmoothingFraction*(filteredAngles-angles);
            smoothSlope = -tan(angles);
        end
        slopeCache = pchip(data(:,1),smoothSlope);
        [breaks,coefs,~,order,dimension] = unmkpp(slopeCache);
        widths = diff(breaks(:));
        integralCoefs = bsxfun(@rdivide,coefs,order:-1:1);
        segmentDrop = sum(integralCoefs.*bsxfun(@power,widths,order:-1:1),2);
        scale = (data(end,2)-data(1,2))/sum(segmentDrop);
        slopeCache = mkpp(breaks,coefs*scale,dimension);
        constants = data(1,2)+[0;cumsum(segmentDrop(1:end-1)*scale)];
        heightCache = mkpp(breaks,[integralCoefs*scale constants],dimension);
    end
    curvatureCache = mkpp(breaks, ...
        bsxfun(@times,coefs(:,1:order-1)*scale,order-1:-1:1),dimension);
    smoothingKey = key;
end
heightPP = heightCache;
slopePP = slopeCache;
curvaturePP = curvatureCache;
end


function [xVisible,yVisible] = trimComparisonSamples(x,y)


firstVisible = find(x >= 0,1,'first');
if isempty(firstVisible)


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


function o = newObstacle()
    o = struct('type', 'Ramp', 'name', '', 'start_m', 20, 'height_m', 0.3, ...
               'up_m', 2, 'flat_m', 0, 'down_m', 2, 'surface', 'Custom', 'mu', 0.04, 'length_m', 5);
end

function o = emptyObstacles()
    o = repmat(newObstacle(), 0, 1);
end

function finish = obstacleEnd(o)
    if strcmp(o.type, 'Ramp')
        span = o.up_m + o.flat_m + o.down_m;
    else
        span = o.length_m;
    end
    finish = o.start_m + span;
end

function validateObstacles(obstacles, courseLength)
    assert(isstruct(obstacles), 'Obstacles must be a struct array.');
    required = fieldnames(newObstacle());
    for k = 1:numel(obstacles)
        o = obstacles(k);
        assert(all(isfield(o, required)), 'An obstacle is missing required fields.');
        assert(ischar(o.type) && any(strcmp(o.type, {'Ramp', 'Surface'})), ...
               'Obstacle type must be Ramp or Surface.');
        assert(ischar(o.name) && (isempty(o.name) || isrow(o.name)), 'Obstacle name must be text.');
        validateattributes(o.start_m, {'numeric'}, {'real', 'finite', 'scalar', 'nonnegative'});
        if strcmp(o.type, 'Ramp')
            validateattributes(o.height_m, {'numeric'}, {'real', 'finite', 'scalar', 'positive'});
            validateattributes(o.up_m, {'numeric'}, {'real', 'finite', 'scalar', 'positive'});
            validateattributes(o.flat_m, {'numeric'}, {'real', 'finite', 'scalar', 'nonnegative'});
            validateattributes(o.down_m, {'numeric'}, {'real', 'finite', 'scalar', 'positive'});
        else
            validateattributes(o.length_m, {'numeric'}, {'real', 'finite', 'scalar', 'positive'});
            validateattributes(o.mu, {'numeric'}, {'real', 'finite', 'scalar', '>=', 0, '<=', 1});
        end
        assert(obstacleEnd(o) <= courseLength, ...
               sprintf('Obstacle %d ends at %.2f m; the course ends at %.2f m.', k, obstacleEnd(o), courseLength));
        for j = 1:k - 1
            sameType = strcmp(o.type, obstacles(j).type);
            overlap = min(obstacleEnd(o), obstacleEnd(obstacles(j))) - max(o.start_m, obstacles(j).start_m);
            assert(~sameType || overlap <= 0, ...
                   sprintf('Obstacles %d and %d overlap. Sections of the same type must not overlap.', j, k));
        end
    end
end


function place(component, row, column)
    component.Layout.Row = row;
    component.Layout.Column = column;
end

function label = addLabel(parent, row, text, bold)
    label = uilabel(parent, 'Text', text, 'FontSize', 12);
    place(label, row, [1 2]);
    if bold
        label.FontWeight = 'bold';
    end
end

function slider = addSlider(parent, row, limits, value, ticks)
    slider = uislider(parent, 'Limits', limits, 'Value', value, 'MajorTicks', ticks);
    place(slider, row, [1 2]);
end

function field = numericControl(parent, row, label, value, limits)
    l = addLabel(parent, row, label, false);
    l.Layout.Column = 1;
    field = uieditfield(parent, 'numeric', 'Limits', limits, 'Value', value);
    place(field, row, 2);
end

function s = onOff(value)
    if value
        s = 'on';
    else
        s = 'off';
    end
end


function validateEditorDraft(o,courseLength)
assert(isstruct(o) && isscalar(o) && all(isfield(o,fieldnames(newObstacle()))), ...
    'Invalid obstacle editor draft.');
assert(ischar(o.type) && any(strcmp(o.type,{'Ramp','Surface'})), 'Invalid draft type.');
assert(ischar(o.name) && (isempty(o.name)||isrow(o.name)), 'Invalid draft name.');
assert(ischar(o.surface) && any(strcmp(o.surface, ...
    {'Custom','Asphalt','Concrete','Grass','Gravel','Sand','Wood'})), 'Invalid surface label.');
keys = {'start_m','height_m','up_m','flat_m','down_m','mu','length_m'};
upper = [courseLength,20,courseLength,courseLength,courseLength,1,courseLength];
for k=1:numel(keys)
    validateattributes(o.(keys{k}),{'numeric'},{'scalar','real','finite','>=',0,'<=',upper(k)});
end
end

function breaks = obstacleBreaks(p)
breaks = [p.xStart p.sectionBoundary1_m p.sectionBoundary2_m p.xEnd];
if isfield(p,'pushEndX_m'), breaks=[breaks p.pushEndX_m]; end
for k=1:numel(p.obstacles)
    o=p.obstacles(k);
    breaks=[breaks o.start_m obstacleEnd(o)];
    if strcmp(o.type,'Ramp')
        breaks=[breaks o.start_m+o.up_m o.start_m+o.up_m+o.flat_m];
    end
end
breaks=unique(breaks(breaks>=p.xStart & breaks<=p.xEnd));
end

function xEnd = findPushEnd(p)


xEnd=p.xStart;
if p.pushForce_N<=0 || p.pushDistance_m<=0, return; end
breaks=obstacleBreaks(p);
remaining=p.pushDistance_m;
for j=1:numel(breaks)-1
    left=breaks(j); right=min(breaks(j+1),left+remaining);
    geometryParameters=p;
    geometryParameters.segmentMid=(breaks(j)+breaks(j+1))/2;
    q=@(x) sqrt(1+hillSlope(x,geometryParameters).^2);
    span=integral(q,left,right,'AbsTol',1e-10,'RelTol',1e-9);
    if abs(remaining-span)<=1e-10
        xEnd=right; return
    elseif remaining<span
        xEnd=fzero(@(x) integral(q,left,x,'AbsTol',1e-10,'RelTol',1e-9)-remaining, ...
            [left right]);
        return
    end
    remaining=remaining-span;
end
xEnd=p.xEnd;
end

function force = startPushForce(x,p)
probe=x;
if isfield(p,'segmentMid'), probe=p.segmentMid+zeros(size(x)); end
force=p.pushForce_N*double(probe>=p.xStart & probe<p.pushEndX_m);
end

function h = rampHeight(x,p)
h=zeros(size(x));
for k=1:numel(p.obstacles)
    o=p.obstacles(k);
    if ~strcmp(o.type,'Ramp'), continue; end
    u=x-o.start_m;
    h=h+o.height_m*max(0,min(min(u/o.up_m,1), ...
        (o.up_m+o.flat_m+o.down_m-u)/o.down_m));
end
end

function dh = rampSlope(x,p)
probe=x;
if isfield(p,'segmentMid'), probe=p.segmentMid+zeros(size(x)); end
dh=zeros(size(x));
for k=1:numel(p.obstacles)
    o=p.obstacles(k);
    if ~strcmp(o.type,'Ramp'), continue; end
    u=probe-o.start_m;
    up=u>=0 & u<o.up_m;
    down=u>=o.up_m+o.flat_m & u<o.up_m+o.flat_m+o.down_m;
    dh(up)=dh(up)+o.height_m/o.up_m;
    dh(down)=dh(down)-o.height_m/o.down_m;
end
end

function mu = surfaceMu(x,p)
probe=x;
if isfield(p,'segmentMid'), probe=p.segmentMid+zeros(size(x)); end
[mu,~]=sectionCoefficients(x,p);
for k=1:numel(p.obstacles)
    o=p.obstacles(k);
    if ~strcmp(o.type,'Surface'), continue; end
    inside=probe>=o.start_m & probe<obstacleEnd(o);
    mu(inside)=o.mu;
end
end

function [mu,b,section] = sectionCoefficients(x,p)
probe=x;


if isfield(p,'segmentMid'), probe=p.segmentMid+zeros(size(x)); end
section=ones(size(x));
section(probe>=p.sectionBoundary1_m)=2;
section(probe>=p.sectionBoundary2_m)=3;
mu=p.Crr1+zeros(size(x)); b=p.b1+zeros(size(x));
mu(section==2)=p.Crr2; b(section==2)=p.b2;
mu(section==3)=p.Crr3; b(section==3)=p.b3;
end

function work = coulombWork(x,p)


starts=[p.xStart p.sectionBoundary1_m p.sectionBoundary2_m];
ends=[p.sectionBoundary1_m p.sectionBoundary2_m p.xEnd];
baseMu=[p.Crr1 p.Crr2 p.Crr3];
weightedDistance=zeros(size(x));
for j=1:3
    left=max(p.xStart,starts(j)); right=min(p.xEnd,ends(j));
    overlap=max(0,min(x,right)-left);
    weightedDistance=weightedDistance+baseMu(j)*overlap;
    for k=1:numel(p.obstacles)
        o=p.obstacles(k);
        if ~strcmp(o.type,'Surface'), continue; end
        surfaceLeft=max(left,o.start_m);
        surfaceRight=min(right,obstacleEnd(o));
        surfaceOverlap=max(0,min(x,surfaceRight)-surfaceLeft);
        weightedDistance=weightedDistance+(o.mu-baseMu(j))*surfaceOverlap;
    end
end
work=p.mass*9.81*weightedDistance;
end

function markObstacles(ax,p,showLabels)
sectionEdges=[p.sectionBoundary1_m p.sectionBoundary2_m];
for j=1:2
    sectionLabel='';
    if showLabels, sectionLabel=sprintf('Section %d',j+1); end
    xline(ax,sectionEdges(j),'--',sectionLabel,'Color',[.45 .49 .54], ...
        'FontSize',9,'HandleVisibility','off');
end
for k=1:numel(p.obstacles)
    o=p.obstacles(k);
    color=[.53 .30 .72];
    label=o.name;
    if strcmp(o.type,'Surface')
        color=[.10 .58 .43];
        label=sprintf('%s (mu=%.3g)',o.name,o.mu);
    end
    if ~showLabels, label=''; end
    xline(ax,o.start_m,':',label,'Color',color,'FontSize',9, ...
        'HandleVisibility','off','Interpreter','none');
    xline(ax,obstacleEnd(o),':','Color',color,'HandleVisibility','off');
end
end


function gps = embeddedGPSSpeed(xEnd)


data = [
    0 0 0;
    0.84221520252094884 0 6.0090000000000003;
    1.6844304050418977 0 12.019;
    2.5266456075628465 0 18.024999999999999;
    16.074325531959474 4.1900000000000004 23.209;
    19.621955236188487 4.1900000000000004 23.295000000000002;
    26.717216090851259 4.1900000000000004 24.305;
    32.238226289982833 4.3499999999999996 25.305;
    37.000468397451975 4.4400000000000004 26.305;
    41.762711627016898 4.25 27.305;
    46.524955976408506 4.25 28.286000000000001;
    50.072590975331615 4.25 29.286000000000001;
    52.090944735656421 3.9900000000000002 30.286000000000001;
    58.828683587875098 3.9399999999999999 31.286000000000001;
    62.376318106500477 3.2799999999999998 32.286000000000001;
    66.413024568363909 3.2799999999999998 33.286000000000001;
    69.173530885509081 2.8900000000000001 34.286000000000001;
    71.934036856107426 2.8900000000000001 35.286000000000001;
    74.694542477090764 2.6200000000000001 36.286000000000001;
    77.221192266222701 2.6200000000000001 37.286000000000001;
    80.76882389571567 2.7599999999999998 38.286000000000001;
    83.529328822087479 2.7599999999999998 39.286000000000001;
    86.055977850934823 2.6899999999999999 40.286000000000001;
    88.816482427817419 2.6899999999999999 41.286000000000001;
    92.182440858314934 2.8199999999999998 42.286000000000001;
    95.730070082466924 3.04 43.286000000000001;
    99.277698824285125 3.1099999999999999 44.286000000000001;
    102.64656166179778 3.1099999999999999 45.286000000000001;
    106.19418992229957 3.3100000000000001 46.284999999999997;
    110.40526783561903 3.54 47.284999999999997;
    114.76067891684805 4.1699999999999999 48.284999999999997;
    119.1160893863568 4.1699999999999999 49.284999999999997;
    123.47149924365094 4.1699999999999999 50.284999999999997;
    128.64568318038681 4.7800000000000002 51.284999999999997;
    133.00109181234845 4.7800000000000002 52.284999999999997;
    135.01944134226881 3.3300000000000001 53.284999999999997;
    138.75644314682663 3.3300000000000001 54.284999999999997;
    141.54625039814053 2.544 55.284999999999997;
    143.77014893133713 2.077 56.284999999999997;
    144.88209819793542 1.552 57.284999999999997;
    145.72431200766837 0.70599999999999996 58.284999999999997;
    146.56652581632986 0.313 59.284999999999997;
    147.96142909644999 0.184 60.284999999999997;
    149.07337836234097 0.377 61.284999999999997;
];
shiftedX = data(:,1)-19.07;
keep = shiftedX > 0 & shiftedX < xEnd;
first = find(keep,1,'first');


endSpeed = interp1(shiftedX,data(:,2),xEnd,'linear',NaN);
endTime = interp1(shiftedX,data(:,3),xEnd,'linear',NaN);
xDisplay = [0;0;shiftedX(keep);xEnd];
vDisplay = [0;data(first,2);data(keep,2);endSpeed];
tDisplay = [NaN;NaN;data(keep,3);endTime];
gps = table(xDisplay,vDisplay,tDisplay, ...
    'VariableNames',{'x_m','speed_mps','time_s'});
end

function speed = sampleRecordedSpeed(gps,x)

[uniqueX,index] = unique(gps.x_m,'last');
speed = interp1(uniqueX,gps.speed_mps(index),x,'linear',NaN);
end
