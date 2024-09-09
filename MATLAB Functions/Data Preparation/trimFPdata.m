function [data] = trimFPdata(data,trimstart,trimend,whichstreams,varargin)
% TRIMDATA    Trims all specified data streams from the index in trimstart to the
%             index in trimend, and adjusts TTLs by the amount trimmed by trimstart.
%
% INPUTS:
%       DATA:           A data frame containing at least the specified
%                       input fields.
%
%       TRIMSTART:      The location to start trimming at.
%
%       TRIMEND:        The location to end trimming at.
%
%       WHICHSTREAMS:   A cell array containing the names of all the
%                       streams to be trimmed.
%  
% OPTIONAL INPUTS:
%       WHICHEPOCS:     A cell array containing the names of all the epocs
%                       to be adjusted due to trimming - subtract the (start
%                       loc - 1) from each specified epoc.
%
% OUTPUTS:
%       DATA:           The data structure with the specified data stream 
%                       containing the trimmed data.
% Written by R M Donka, February 2024
% Stored in Roitman Photometry GitHub repository, see Wiki for additional notes.
    inputs = struct(...
        'whichepocs',[]);
    inputs = parseArgsLite(varargin,inputs);

    if isempty(inputs.whichepocs) == true
        disp('TRIMMING DATA: No epocs specified for adjustment.')
    end

    for eachfile = 1:length(data)
        disp(append('   TRIMMING FILE NUMBER ',num2str(eachfile)))
    
        startloc = data(eachfile).(trimstart); % start loc
        endloc = data(eachfile).(trimend);
    
        for eachstream = 1:length(whichstreams)
            stream = char(whichstreams(eachstream));
            try 
                data(eachfile).(stream) = data(eachfile).(stream)(startloc:endloc);
            catch
                disp(append('WARNING: File number ',num2str(eachfile), ' - failed to trim stream: ', stream))
            end
        end
        
        if isempty(inputs.whichepocs) == false
            for eachepoc = 1:length(whichepocs)
                epoc = char(whichepocs(eachepoc));
                try 
                    data(eachfile).(epoc) = data(eachfile).(epoc)-(startloc-1);
                catch
                    disp(append('WARNING: File number ',num2str(eachfile), ' - failed to adjust epoc: ', epoc))
                end
            end
        end
    end
end