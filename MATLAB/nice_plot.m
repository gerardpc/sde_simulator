%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nice_plot(xx, yy, x_label, y_label, t_title, varargin)
%
% Plot vectors nicely, in "journal" fashion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%           xx       -- x vector
%           yy       -- y vector
%           x_label  -- x axis label
%           y_label  -- y axis label
%           t_title  -- Plot title
%           varargin -- color, linestyle, linewidth, optional inputs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GP Conangla, 7.2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[] = nice_plot(xx, yy, x_label, y_label, t_title, varargin)

%% Varargins
% color specified?
specify_color = find(ismember(varargin, 'color'));
if(~isempty(specify_color))
    color = varargin{specify_color + 1};
else
    Fig_num = get(gcf, 'Number');
    current_count = get(Fig_num, 'UserData');
    if isempty(current_count)
        current_count = 6;
    end
    color_list = get(gca,'ColorOrder');
    color = color_list(mod(current_count, 7) + 1, :);
    current_count = current_count + 1;
    set(Fig_num, 'UserData', current_count);
end

% linestyle specified?
specify_linestyle = find(ismember(varargin, 'linestyle'));
if(~isempty(specify_linestyle))
    linestyle = varargin{specify_linestyle + 1};
else
    linestyle = '-';
end

% linewidth specified?
specify_linewidth = find(ismember(varargin, 'linewidth'));
if(~isempty(specify_linewidth))
    linewidth = str2double(varargin{specify_linewidth + 1});
else
    linewidth = 2;
end

% linewidth specified?
specify_size = find(ismember(varargin, 'size'));
if(~isempty(specify_size))
    marker_size = str2double(varargin{specify_size + 1});
else
    marker_size = 40;
end

%% Now plot
% plot or scatter?
plot_type = find(ismember(varargin, 'scatter'));
if(~isempty(plot_type))    
    % Edge or not edge?
    plot_type = find(ismember(varargin, 'edge'));    
    if(~isempty(plot_type))
        scatter(xx, yy, marker_size, 'MarkerFaceColor', color, 'MarkerFaceAlpha', 1, 'MarkerEdgeColor', 'black');
        hold on;  
    else         
        scatter(xx, yy, marker_size, 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
        hold on;  
        area(xx, yy, 'EdgeAlpha', 0, 'FaceColor', color, 'FaceAlpha', 0.2);
    end
    
    % linestyle specified?
    plot_type = find(ismember(varargin, 'semilogy'));
    if(~isempty(plot_type))
        set(gca, 'YScale', 'log');
    end
    plot_type = find(ismember(varargin, 'semilogx'));
    if(~isempty(plot_type))
        set(gca, 'XScale', 'log');
    end
    plot_type = find(ismember(varargin, 'loglog'));
    if(~isempty(plot_type))
        set(gca, 'XScale', 'log');
        set(gca, 'YScale', 'log');
    end
else
    plot(xx, yy, 'LineStyle', linestyle, 'color', color, 'LineWidth', linewidth);
    hold on;
    plot_type = find(ismember(varargin, 'semilogy'));
    if(~isempty(plot_type))
        set(gca, 'YScale', 'log');
    end
    plot_type = find(ismember(varargin, 'semilogx'));
    if(~isempty(plot_type))
        set(gca, 'XScale', 'log');
    end
    plot_type = find(ismember(varargin, 'loglog'));
    if(~isempty(plot_type))
        set(gca, 'XScale', 'log');
        set(gca, 'YScale', 'log');
    end
end

%% Nice labels, titles and sizes
set(gca,'FontSize',18,'FontName', 'CMU Sans Serif');
xlabel(x_label, 'fontsize', 24, 'FontName', 'CMU Sans Serif', 'interpreter', 'latex');
ylabel(y_label, 'fontsize', 24, 'FontName', 'CMU Sans Serif', 'interpreter', 'latex');
t_title = strrep(t_title, '_', ' '); % clean title string '_' -> ' ' 
title(t_title);
xlim([min(xx), max(xx)]);
grid on;
ax = gca;
ax.GridAlpha = 0.3;
box on;

end