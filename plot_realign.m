function P = plot_realign(P)
% Plot reaigment parameters
% FORMAT P = plot_realign(P)
%
% P     - char array of filenames
%         All operations are performed relative to the first image.
%         ie. Coregistration is to the first image, and resampling
%         of images is into the space of the first image.
%         For multiple sessions, P should be a cell array, where each
%         cell should be a matrix of filenames.
%-Images
%--------------------------------------------------------------------------
if ~iscell(P), P = {P}; end
for i=1:numel(P), if ischar(P{i}), P{i} = spm_vol(P{i}); end; end
P(cellfun(@isempty,P)) = [];

P = cat(1,P{:});
if length(P)<2, return; end
Params = zeros(numel(P),12);
for i=1:numel(P)
    Params(i,:) = spm_imatrix(P(i).mat/P(1).mat);
end

%-Display results: translation and rotation over time series
%--------------------------------------------------------------------------
fg = figure(101);
clf;
ax = axes('Position',[0.1 0.65 0.8 0.2],'Parent',fg,'Visible','off');
set(get(ax,'Title'),'String','Image realignment',...
    'FontSize',16,'FontWeight','Bold','Visible','on');
x     =  0.1;
y     =  0.9;
for i = 1:min([numel(P) 4])
    text(x,y,[sprintf('%-4.0f',i) P(i).fname ',' num2str(P(i).n(1))],...
        'FontSize',10,'Interpreter','none','Parent',ax);
    y = y - 0.08;
end
if numel(P) > 12
    text(x,y,'................ etc','FontSize',10,'Parent',ax); end

ax = axes('Position',[0.1 0.35 0.8 0.2],'Parent',fg,'XGrid','on','YGrid','on',...
    'NextPlot','replacechildren','ColorOrder',[0 0 1;0 0.5 0;1 0 0]);
plot(Params(:,1:3),'Parent',ax)
s  = {'x translation','y translation','z translation'};
%text([2 2 2], Params(2, 1:3), s, 'Fontsize',10,'Parent',ax)
legend(ax, s, 'Location','Best')
set(get(ax,'Title'),'String','translation','FontSize',16,'FontWeight','Bold');
set(get(ax,'Xlabel'),'String','image');
set(get(ax,'Ylabel'),'String','mm');


ax = axes('Position',[0.1 0.05 0.8 0.2],'Parent',fg,'XGrid','on','YGrid','on',...
    'NextPlot','replacechildren','ColorOrder',[0 0 1;0 0.5 0;1 0 0]);
plot(Params(:,4:6)*180/pi,'Parent',ax)
s  = {'pitch','roll','yaw'};
%text([2 2 2], Params(2, 4:6)*180/pi, s, 'Fontsize',10,'Parent',ax)
legend(ax, s, 'Location','Best')
set(get(ax,'Title'),'String','rotation','FontSize',16,'FontWeight','Bold');
set(get(ax,'Xlabel'),'String','image');
set(get(ax,'Ylabel'),'String','degrees');


