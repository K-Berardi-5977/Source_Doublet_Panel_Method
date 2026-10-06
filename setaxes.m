function setaxes
% SETAXES  Apply standardized axis styling to the current axes.
% No inputs required. Operates on gca.

ax = gca;

% --- Font settings ---
scale = 1.5;
baseFont = ax.FontSize;

ax.FontName = 'Times New Roman';
ax.FontSize = baseFont * scale;

% --- Grid settings ---
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = '--';

% Major grid (bold)
ax.GridAlpha = 0.35;
ax.LineWidth = 1.0;

% Minor grid (lighter)
ax.MinorGridAlpha = 0.15;

% --- Title and axis labels ---
t = get(ax,'Title');
xl = get(ax,'XLabel');
yl = get(ax,'YLabel');

set([t xl yl], ...
    'FontWeight','bold', ...
    'FontName','Times New Roman', ...
    'FontSize', baseFont * scale);
end