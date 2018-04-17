function formatPlot(fig, varargin)
% formats plot to a publishable format.
% use with export_fig
% 
% Options:
%	'-zbuffer'        : use zbuffer as renderer
%	'twocolumn-single': optimized for 8.6cm x 5.36cm resizing
%
% by Florian Schroeder 04/03/2014
%   Modified:
%       - FS 01/09/2014: 'zbuffer' uses smaller font sizes
%       - FS 10/09/2014: restricted findobj() to gca! (much faster)

%% Preferences
f = figure(fig);        % activate figure
subplots = get(gcf, 'children');    % array of subplot references.

% change size
screen = get(0,'screensize');	% x,y,width,height
if f.WindowStyle == 'normal'
	f.Position = [min(screen([3,4])-[1280,800]-100,f.Position([1,2])),1280 800];	% 100 for title bar
	% creates for 40pt: 31.5cm x 19.63cm @96dpi
	% creates for 21pt: 29.1cm x 19.31cm @96dpi
end

plFontSize = 40;
if ismember(varargin,'zbuffer')
    plFontSize = 35;
end
plLineWidth = 4;	%axes have width!
if ismember(varargin,'twocolumn-single')		% seems to work!
	% intended to be resized to 8.6cm x 5.36cm: 
	%     for 31.5cm -> 350dpi. -> factor 3.6628
	%     for 29.1cm -> 325dpi. -> factor 3.3854
	% needs '-m2' for 700dpi in exportfig
	% want 2mm Font Height on paper.
	% Font measured in 72dpi = 28.35dpcm -> FontSize[pt] = 28.35 * 3.6628 * height[cm]
	plFontSize = 20;					% = 2.022mm
	plLineWidth = 4;					% = 0.385mm
end
plAxesLineWidth = plLineWidth/2;
colors=lines(3);

for i = 1:1:length(subplots)
    gca = subplots(i);              % activate subplot
    
    % get(get(gcf,'Children'),'Type') to see all types.
    if ~strcmp(get(gca,'Type'),'axes')
        continue;
    end
    
    % get title and reset fontsize!
    titleref = get(gca,'title');
    set(titleref,'FontSize',plFontSize);

    xref = get(gca,'xlabel');
    set(xref,'FontSize',plFontSize);

    yref = get(gca,'ylabel');
    set(yref,'FontSize',plFontSize);

    zref = get(gca,'zlabel');
    set(zref,'FontSize',plFontSize);

    set(gca,'LineWidth',plAxesLineWidth, 'FontSize',plFontSize);
    box(gca,'on');

    h=findobj(gca,'Type', 'line');
    set(h,'LineWidth',plLineWidth);
    
    h=findobj(gca,'Type', 'text');
    set(h,'FontSize',plFontSize);

end

end