function [data] = trimFPdata(data,whichtrimstart,whichtrimend,whichstreams,varargin)
% TRIMDATA    Trims all specified data streams from the index in whichtrimstart to the
%             index in whichtrimend, and adjusts TTLs by the amount trimmed by whichtrimstart.
%
% INPUTS:
%       DATA:           Data structure; A data structure containing at 
%                       least the specified input fields.
%
%       WHICHTRIMSTART: String; The name of the field containing the 
%                       location of the start of the session; everything
%                       before will be trimmed. For example: 'sessionstart'
%
%       WHICHTRIMEND:   String; The name of the field containing the 
%                       location of the end of the session; everything
%                       after will be trimmed. For example: 'sessionend'
%
%       WHICHSTREAMS:   Cell array; Contains the names of all the
%                       streams to be trimmed.
%                       For example: whichstreams = {'sig', 'baq'};
% OPTIONAL INPUTS:
%       WHICHEPOCS:     A cell array containing the names of all the epocs
%                       to be adjusted due to trimming - subtract the (start
%                       loc - 1) from each specified epoc.
%                       For example: whichepocs = {'injt', 'trialstart', 'trialend'};
%
% OUTPUTS:
%       DATA:           Data structure; The original data structure with 
%                       the specified data streams containing the trimmed 
%                       data and the specified epocs adjusted.
%
% Written by R M Donka, February 2024
% Stored in the PASTa GitHub Repository, see the user guide for additional
% documentation: https://rdonka.github.io/PASTa/

    
    %% Prepare Settings
    % Import required and optional inputs into a structure
    inputs = struct(...
        'whichepocs',[]);

    myinputs = parseArgsLite(varargin,inputs);

    if isempty(inputs.whichepocs) == true
        disp('TRIMMING DATA: No epocs specified for adjustment.')
    else
        disp('TRIMMING DATA: Specified epocs will be adjusted.')
        disp('EPOCS:') % Display specified epocs
        disp(inputs.whichepocs)
    end

    %% Trim data
    for eachfile = 1:length(data)
        disp(append('   TRIMMING FILE NUMBER ',num2str(eachfile)))
    
        startloc = data(eachfile).(whichtrimstart); % Extract session start location - everything before will be trimmed
        endloc = data(eachfile).(whichtrimend); % Extract session end location - everything after will be trimmed
    
        for eachstream = 1:length(whichstreams) % Trim each stream
            stream = char(whichstreams(eachstream));
            try 
                data(eachfile).(stream) = data(eachfile).(stream)(startloc:endloc);
            catch
                disp(append('WARNING: File number ',num2str(eachfile), ' - failed to trim stream: ', stream))
            end
        end
        
        if isempty(inputs.whichepocs) == false % If WHICHEPOCS is used, adjust each epoc by the number of samples trimmed at the start of the session
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