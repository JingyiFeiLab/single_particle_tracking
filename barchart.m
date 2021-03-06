% Plots a bar chart and gives a different color to each bar.
% Also plots the value of the bar above the bar.
clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables.
workspace;  % Make sure the workspace panel is showing.
fontSize = 30;
format compact

% Ask user for the number of bars.
defaultValue = 7;
titleBar = 'Enter an integer value';
userPrompt = 'Enter the number of bars';
caUserInput = inputdlg(userPrompt, titleBar, 1, {num2str(defaultValue)});
if isempty(caUserInput),return,end; % Bail out if they clicked Cancel.
integerValue = round(str2double(cell2mat(caUserInput)));
% Check for a valid integer.
if isnan(integerValue)
    % They didn't enter a number.  
    % They clicked Cancel, or entered a character, symbols, or something else not allowed.
    integerValue = defaultValue;
    message = sprintf('I said it had to be an integer.\nI will use %d and continue.', integerValue);
    uiwait(warndlg(message));
end

% Define sample data in the range 20-80.
x = 1 : integerValue;
y = 20 + 80 * rand(integerValue)
numberOfBars = length(y);

button = menu('Use which colormap?', 'Custom', 'Random', 'Jet', 'Hot', 'Lines');
if button == 1
	% Make up a custom colormap specifying the color for each bar series.
	barColorMap(1,:) = [.2 .71 .3];	% Green Color for segment 1.
	barColorMap(2,:) = [.25 .55 .79];	% Blue Color for segment 2.
	barColorMap(3,:) = [.9 .1 .14];	% Red Color for segment 3.
	barColorMap(4,:) = [.9 .9 .14];	% Yellow Color for segment 4.
	% I have not defined any more than 4 colors in this demo.
	% For any number of bars beyond 4, just make up random colors.
	if numberOfBars > 4
		barColorMap(5:numberOfBars, 1:3) = rand(numberOfBars-4, 3);
	end
elseif button == 2
	% Example of using colormap with random colors
	barColorMap = rand(numberOfBars, 3); 
elseif button == 3
	% Example of using pre-defined jet colormap
	barColorMap = jet(numberOfBars); 
elseif button == 4
	% Example of using pre-defined Hot colormap
	barColorMap = hot(numberOfBars); 
else
	% Example of using pre-defined lines colormap
	barColorMap = lines(numberOfBars); 
end

% Plot each number one at a time, calling bar() for each y value.
barFontSize = 15;
for b = 1 : numberOfBars
	% Plot one single bar as a separate bar series.
	handleToThisBarSeries(b) = bar(x(b), y(b), 'BarWidth', 0.9);
	% Apply the color to this bar series.
	set(handleToThisBarSeries(b), 'FaceColor', barColorMap(b,:));
	% Place text atop the bar
	barTopper = sprintf('y(%d) = %.3f', x(b), y(b));
	text(x(b)-0.2, y(b)+3, barTopper, 'FontSize', barFontSize);
	hold on;
end

% Fancy up the graph.
grid on;
caption = sprintf('Data plotted in %d barseries, each with a different color', length(y));
title(caption, 'FontSize', fontSize);
xlabel('x', 'FontSize', fontSize);
ylabel('y', 'FontSize', fontSize);
% Restore the x tick marks.
set(gca, 'XTickMode', 'Auto');
% set(gca, 'XTickLabels', xTickLabels);
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% Give a name to the title bar.
set(gcf,'name','Demo by ImageAnalyst','numbertitle','off');


