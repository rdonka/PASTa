function [alltraces] = plotTraces(data,whichfile,whichstreams,whichtitles,whichylabels,whichcolors)
% PLOTTRACES    Plots whole session fiber photometry traces. This function
% will plot every stream input in whichstreams. Make sure input cell arrays
% are the same length (ie, each stream has a title, y lable, and color.
%
% INPUTS:
%       DATA:           This is a structure that contains at least the
%                       streams to be plotted.
%
%       WHICHFILE:      The file number to plot (this can be set in a for
%                       loop to plot all files).
%
%       WHICHSTREAMS:   A cell array containing the field names of the
%                       streams to be plotted. Reccomended to include raw, 
%                       subtracted, subtracted and filtered, and normalized 
%                       data.
%
%       WHICHTITLES:    A cell array containing the titles to be added to
%                       each stream sub plot.
%
%       WHICHYLABELS:   A cell array containing the y axis labels to be
%                       added to each stream sub plot.
%
%       WHICHCOLORS:    A cell array containing the html color codes to be
%                       used for each stream sub plot.
%
% OUTPUT:
%       ALLTRACES:      A plot object containing subplots for each input
%                       stream.
% Written by R M Donka, March 2024.
% Stored in RoitmanPhotometry GitHub repository, see Wiki for additional notes.
    close all
    ntraces = length(whichstreams);

    alltraces = figure();
    tiledlayout(ntraces, 1, 'Padding','compact', 'TileSpacing','compact');
    alltraces.Units = 'inches';
    alltraces.OuterPosition = [0.25 0.25 5 ntraces*1.5];
    
        for eachplot = 1:length(whichstreams)
            currstream = char(whichstreams(eachplot));
            currxlength = length(data(whichfile).(currstream));
            currtitle = char(whichtitles(eachplot));
            currYlabel = char(whichylabels(eachplot));
            currcolor = char(whichcolors(eachplot));
        
            alltraces(eachplot) = nexttile;
            hold on
            plot(data(whichfile).(currstream),'Color',currcolor)
            xlim([0 currxlength])
            title(currtitle)
            xlabel('Sample')
            ylabel(currYlabel)
            hold off
            
        end
    
    fontsize(gcf,scale=.8)
end