function [data] = findsubjectmeantrialtrace(data,trialstreamfieldname,trialtypefieldname,trialeventepocfieldnames,varargin)
% FINDSUBJECTMEANTRIALTRACE    Finds the mean trace by trial type for each subject.
%
% INPUTS:
%   DATA:               This is a structure that contains at least the
%                           specified field with trial by trial signal and 
%                           trial type.
%
%   TRIALSTREAMFIELDNAME:   The name of the field in the data structure that
%                           contains the cut trial signal.
%       
%   TRIALTYPEFIELDNAME:  The name of the field in the data structure
%                            that contains the trial type identifier (trials of the same type
%                            will be averaged).
%
%   TRIALEVENTEPOCFIELDNAMES: Cell array; names of fields containing
%                               relevant trial epochs. These will be averaged and rounded to
%                           identify common mean event epochs.
%
% OPTIONAL INPUTS:
%   TRIALIDOUPUTFIELDNAME: String; Output name for the trial ID field
%                          output with the mean trace data. Default:
%                          "TrialType"
% OUTPUTS:
%       DATA:           This is the original data structure with added 
%                       fields for the mean signal for each trial type.
%                       Fields will be named as 
%                       data.mean_TRIALSTREAMFIELDNAME with each trial type in the field "TrialType".
% 
% Written by R M Donka, January 2024
% Stored in RoitmanPhotometry GitHub repository, see Wiki for additional notes.

    %% Prepare Settings
    % Import required and optional inputs into a structure
    p = createParser(mfilename); % Create parser object with custom settings - see createParser helper function for more details
    addParameter(p, 'trialidoutputfieldname', 'TrialType'); 

    parse(p, varargin{:});

    % Retrieve parsed inputs into params structure
    params = p.Results;

    % Initialize params
    params.trialstreamfieldname = trialstreamfieldname;
    params.trialtypefieldname = trialtypefieldname;
    params.trialeventepocfieldnames = trialeventepocfieldnames;
    
    if ~isempty(params.trialidoutputfieldname)
        trialidoutputfieldname = params.trialidoutputfieldname;
    end

    outputtrialstreamfieldname = append('mean_',trialstreamfieldname);
    for eachfile = 1:length(data)
        fprintf('Find trial means for file number: %.f \n',eachfile) % Display which file is being subtracted
        trialtypes = unique([data(eachfile).(trialtypefieldname)]);
        trialmeanstruct = struct();
        try
            for eachtrialtype = 1:length(trialtypes) % Save each trial type to a row
                currtrialtype = trialtypes(eachtrialtype);
                currtrialtypeidxs = find([data(eachfile).(trialtypefieldname)] == currtrialtype);


                trialmeanstruct.(trialidoutputfieldname)(eachtrialtype,1) = trialtypes(eachtrialtype);
                trialmeanstruct.trialdata(eachtrialtype,:) = mean(data(eachfile).(trialstreamfieldname)(currtrialtypeidxs,:));

                for eachtrialeventepoc = 1:length(trialeventepocfieldnames)
                    currepocname = trialeventepocfieldnames{eachtrialeventepoc};
                    trialmeanstruct.(currepocname)(eachtrialtype,1) = mean(data(eachfile).(currepocname)(currtrialtypeidxs,:));
                end
            end

            data(eachfile).(outputtrialstreamfieldname) = trialmeanstruct; % Add to data structure
        catch
            disp(['Find trial means without success:' num2str(eachfile)]);
        end
    end
end