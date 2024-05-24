function [] = PlotOri(Node,Panel,Trigl,lstyle,alfa,color)
if nargin < 4
    lstl = '-'; al = 1; cc = 'g';
else
    lstl = lstyle; al = alfa; cc = color;
end
if ~isempty(Trigl)
    patch('faces', Trigl, 'vertices', Node, 'facecolor', cc, ...
          'linestyle', 'none', 'facelighting', 'flat', 'edgecolor', (1-al)*[1 1 1]);
end
hold on;
Panelsize = cellfun(@numel,Panel);
Ptri = cell(sum(Panelsize==3),1);
Pquad = cell(sum(Panelsize==4),1);
flg3 = find(Panelsize==3);
flg4 = find(Panelsize==4);
for i = 1:numel(flg3), Ptri{i} = Panel{flg3(i)}; end
for j = 1:numel(flg4), Pquad{j} = Panel{flg4(j)}; end

if ~isempty(Trigl)
    patch('faces', cell2mat(Ptri), 'vertices', Node, 'facecolor', 'none', ...
          'linestyle', lstl, 'linewidth', 1, 'edgecolor', (1-al)*[1 1 1]);
    patch('faces', cell2mat(Pquad), 'vertices', Node, 'facecolor', 'none', ...
          'linestyle', lstl, 'linewidth', 1, 'edgecolor', (1-al)*[1 1 1]);
else
    patch('faces', cell2mat(Ptri), 'vertices', Node, 'facecolor', cc, ...
          'linestyle', lstl, 'linewidth', 1, 'edgecolor', (1-al)*[1 1 1]);
    patch('faces', cell2mat(Pquad), 'vertices', Node, 'facecolor', cc, ...
          'linestyle', lstl, 'linewidth', 1, 'edgecolor', (1-al)*[1 1 1]);
end