function lgd = addCheckboxLegend(ax)
% addCheckboxLegend  Add a checkbox-style legend to an axes
%
%   lgd = addCheckboxLegend(ax)
%
% - Uses Line objects in AX
% - Excludes the last line (assumed baseline)
% - Uses DisplayName for labels
% - Clicking legend entries toggles visibility
%
% Returns legend handle LGD

    if nargin < 1 || isempty(ax)
        ax = gca;
    end

    % Get line objects in creation order (oldest -> newest)
    ch = flipud(ax.Children);
    isLine = arrayfun(@(h) isa(h,'matlab.graphics.chart.primitive.Line'), ch);
    lines = ch(isLine);

    if numel(lines) < 2
        error('Need at least 2 lines (conditions + baseline).');
    end

    % Exclude last line (baseline)
    condLines = lines(1:end-1);

    % Get labels
    names = get(condLines,'DisplayName');
    if ischar(names), names = {names}; end

    for k = 1:numel(names)
        if isempty(names{k})
            names{k} = sprintf('Curve %d',k);
        end
    end

    % Initial visibility state
    isOn = arrayfun(@(h) strcmp(h.Visible,'on'), condLines);

    % Create legend strings with checkboxes
    labels = addCheckboxes(names, isOn);

    % Create legend
    lgd = legend(ax, condLines, labels, ...
        'Interpreter','none', ...
        'Location','eastoutside');

    % Store state
    lgd.UserData.lines = condLines;
    lgd.UserData.names = names;
    lgd.UserData.isOn  = isOn;

    % Click-to-toggle
    lgd.ItemHitFcn = @(src,evt) toggleLine(src, evt);
end

% ---------- helpers ----------
function out = addCheckboxes(baseLabels, isOn)
    out = baseLabels;
    for k = 1:numel(baseLabels)
        if isOn(k)
            out{k} = ['☑ ' baseLabels{k}];
        else
            out{k} = ['☐ ' baseLabels{k}];
        end
    end
end

function toggleLine(lgd, evt)
    L = lgd.UserData.lines;

    idx = find(L == evt.Peer, 1, 'first');
    if isempty(idx), return; end

    newState = ~strcmp(L(idx).Visible,'on');
    L(idx).Visible = ternary(newState,'on','off');

    lgd.UserData.isOn(idx) = newState;
    lgd.String = addCheckboxes(lgd.UserData.names, lgd.UserData.isOn);
end

function out = ternary(cond,a,b)
    if cond, out = a; else, out = b; end
end



