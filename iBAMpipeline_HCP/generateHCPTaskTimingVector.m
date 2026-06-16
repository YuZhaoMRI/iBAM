function TaskTimings = generateHCPTaskTimingVector(events_file, TR, Descent_dur, dumypoints)
%GENERATETASKTIMINGVECTOR Creates a binary task timing vector from an event file.
%
%   TaskTimings = GENERATETASKTIMINGVECTOR(events_file, TR, Descent_dur, dumypoints)
%   generates a vector indicating task-on periods based on event onsets
%
%   Input:
%       events_file  - Path to the event file (plain text). The file should
%                      contain an N×2 matrix: column 1 = onset time (seconds),
%                      column 2 = duration (seconds).
%       TR           - Repetition time of the fMRI scan (seconds).
%       Descent_dur  - Duration of an additional post-event period (seconds),
%    
%       dumypoints   - Number of initial volumes to exclude/shift (e.g., for
%                      removing dummy scans). The task vector starts earlier
%                      by this amount to allow later alignment.
%
%   Output:
%       TaskTimings  - A vector [1 × total_volumes] where -1 indicates
%                      neutual task-on periods; 1 indicates
%                      fear task-on periods;and 0 indicates rest.
%
%
%   Example:
%       TaskTimings = generateTaskTimingVector('events.txt', 0.8, 10, 5);


    % Load event file (assumed to be a plain text file with N×2 matrix)
    events_1 = load([events_file, 'neut.txt']);  % N×2 matrix
    events_2 = load([events_file, 'fear.txt']);  % N×2 matrix
    % Column 1: onset time (seconds), Column 2: duration (seconds)
    
    % Determine total scan duration (end of last event + extra margin)
    t_end = events_2(2,1) + events_2(2,2) + 2*Descent_dur;
    
    % Initialize binary timing vector (one element per volume)
    TaskTimings = zeros(1, ceil(t_end/TR + 0.5));
    t_N = size(TaskTimings, 2);  % Total number of volumes (unused here but kept for clarity)
    
    % Fill neut task-on periods
    for i = 1:2%size(events, 1)
        % Convert onset to volume index and shift left by dummy points
        start_t = floor(events_1(i,1) / TR) - dumypoints;
        
        % Convert duration to number of volumes
        dur = floor(events_1(i,2) / TR);
        
        % Indices for this block
        idx = start_t : (start_t + dur - 1);  % Note: dur volumes means dur elements
        
        % Set task-on (protect against out-of-bounds indices)
        valid_idx = idx(idx >= 1 & idx <= length(TaskTimings));
        TaskTimings(valid_idx) = -1;
    end

        % Fill fear task-on periods
    for i = 1:2%size(events, 1)
        % Convert onset to volume index and shift left by dummy points
        start_t = ceil(events_2(i,1) / TR) - dumypoints;
        
        % Convert duration to number of volumes
        dur = ceil(events_2(i,2) / TR);
        
        % Indices for this block
        idx = start_t : (start_t + dur - 1);  % Note: dur volumes means dur elements
        
        % Set task-on (protect against out-of-bounds indices)
        valid_idx = idx(idx >= 1 & idx <= length(TaskTimings));
        TaskTimings(valid_idx) = 1;
    end

end

