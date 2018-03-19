clear all
close all
clc

h = plot(randn(7));
for i = 1:7
 colors(i,:) = h(i).Color;
end
close all

data = load('pals_data/extracted/Howard.txt');

ET = data(:,9);
H  = data(:,10);
Ta = data(:,5);
Rn = data(:,6);

I = find(isnan(data(:,1)));
ET(I) = [];
H(I) = [];
Ta(I) = [];
Rn(I) = [];

X = 1:length(Ta);

I = 1:300;

figure(1); close(1); figure(1);
set(gcf,'color','w','position',[1,1,1500,1000]);

subplot(2,1,1)
[ax,h1,h2] = plotyy(X(I),Ta(I),X(I),Rn(I));
h1.LineWidth = 4;
h2.LineWidth = 4;
h1.Color = colors(1,:);
h2.Color = colors(2,:);
ax(1).YLabel.FontSize = 16;
ax(2).YLabel.FontSize = 16;
ax(1).YLabel.String = 'Air Temperature [K]';
ax(2).YLabel.String = 'Net Radiation [W/m^2]';
ax(1).XLim = [X(I(1)),X(I(end))];
ax(2).XLim = [X(I(1)),X(I(end))];
ax(1).YColor = colors(1,:);
ax(2).YColor = colors(2,:);
title('Model Inputs','fontsize',20);

subplot(2,1,2)
[ax,h1,h2] = plotyy(X(I),ET(I),X(I),H(I));
h1.LineWidth = 4;
h2.LineWidth = 4;
h1.Color = colors(3,:);
h2.Color = colors(4,:);
ax(1).YLabel.FontSize = 16;
ax(2).YLabel.FontSize = 16;
ax(1).YLabel.String = 'Latent Heat [W/m^2]';
ax(2).YLabel.String = 'Sensible Heat [W/m^2]';
ax(1).XLabel.String = 'Timestep';
ax(1).XLabel.FontSize = 16;
ax(1).XLim = [X(I(1)),X(I(end))];
ax(2).XLim = [X(I(1)),X(I(end))];
ax(1).YColor = colors(3,:);
ax(2).YColor = colors(4,:);
title('Model Outputs','fontsize',20);

