function sample = time2sample(ts, varargin)

% defaults
FS  = 195312.5;
T1  = 0;
T2 = 0;
TO_TIME = 0;

VALID_PARS = {'FS','T1','T2','TO_TIME'};

% parse varargin
for ii = 1:2:length(varargin)
    if ~ismember(upper(varargin{ii}), VALID_PARS)
        error('%s is not a valid parameter. See help time2sample.', upper(varargin{ii}));
    end
    eval([upper(varargin{ii}) '=varargin{ii+1};']);
end

sample = ts * FS;
if T2
    % drop precision beyond 1e-9
    exact = round(sample * 1e9) / 1e9;
    sample = floor(sample);
    if exact == sample
        sample = sample - 1;
    end
else
    % drop precision beyond 1e-9
    sample = round(sample * 1e9) / 1e9;
    if T1
        sample = ceil(sample);
    else
        sample = round(sample);
    end
sample = double(sample);
if TO_TIME
    sample = double(sample) / FS;
end
end