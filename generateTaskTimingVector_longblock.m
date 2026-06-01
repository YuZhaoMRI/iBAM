function [TaskTimings, TaskTimings2]= generateTaskTimingVector_longblock(events_file, TR, Descent_dur, dumypoints)
%GENERATETASKTIMINGVECTOR Creates a binary task timing vector from an event file.
%
%   TaskTimings = GENERATETASKTIMINGVECTOR(events_file, TR, Descent_dur, dumypoints)
%   generates a binary vector indicating task-on periods based on event onsets
%   and durations stored in a text file. This vector is commonly used for
%   task regressor creation in GLM analysis of task-based fMRI data.
%
%   Input:
%       events_file  - Path to the event file (plain text). The file should
%                      contain an N×2 matrix: column 1 = onset time (seconds),
%                      column 2 = duration (seconds).
%       TR           - Repetition time of the fMRI scan (seconds).
%       Descent_dur  - Duration of an additional post-event period (seconds),
%                      often used for hemodynamic response tail modeling.
%       dumypoints   - Number of initial volumes to exclude/shift (e.g., for
%                      removing dummy scans). The task vector starts earlier
%                      by this amount to allow later alignment.
%
%   Output:
%       TaskTimings  - Binary vector [1 × total_volumes] where 1 indicates
%                      task-on periods and 0 indicates rest.
%
%   Notes:
%       - The total length is determined by the end of the last event plus
%         a safety margin (2 × Descent_dur).
%       - Onsets and durations are converted to volume indices using ceil().
%       - The vector is zero-initialized and filled with 1s during task blocks.
%
%   Example:
%       TaskTimings = generateTaskTimingVector('events.txt', 0.8, 10, 5);


    % Load event file (assumed to be a plain text file with N×2 matrix)
    events = load(events_file);  % N×2 matrix
    
    % Column 1: onset time (seconds), Column 2: duration (seconds)
    
    % Determine total scan duration (end of last event + extra margin)
    t_end = max(events(:,1) + events(:,2) + 2*Descent_dur);
    
    % Initialize binary timing vector (one element per volume)
    TaskTimings = zeros(1, ceil(t_end/TR + 0.5));
    TaskTimings2 = zeros(1, ceil(t_end/TR + 0.5));
    t_N = size(TaskTimings, 2);  % Total number of volumes (unused here but kept for clarity)
    
    long_dur=quantile(events(:,2),06);
    % Fill task-on periods
    for i = 1:size(events, 1)
        % Convert onset to volume index and shift left by dummy points 

        if (i==18||i==23||i==25||i==27||i==28)~=1
            continue
        end

        if events(i,2)<long_dur
            continue
        end
        start_t = floor(events(i,1) / TR) - dumypoints;
       
        % Convert duration to number of volumes
        dur = floor(events(i,2) / (2*TR));
        
        % Indices for this block
        idx = start_t : (start_t + dur - 1);  % Note: dur volumes means dur elements
        
        % Set task-on (protect against out-of-bounds indices)
        valid_idx = idx(idx >= 1 & idx <= length(TaskTimings));
        TaskTimings(valid_idx) = 1;
    end


    for i = 1:size(events, 1)
        % Convert onset to volume index and shift left by dummy points 
        if (i==18||i==23||i==25||i==27||i==28)~=1
            continue
        end

        if events(i,2)<long_dur
            continue
        end
        start_t = floor(events(i,1) / TR) - dumypoints;
       
        % Convert duration to number of volumes
        dur = floor(events(i,2) / (2*TR));
        
        % Indices for this block
        idx = start_t+dur: (start_t + 2*dur - 1);  % Note: dur volumes means dur elements
        
        % Set task-on (protect against out-of-bounds indices)
        valid_idx = idx(idx >= 1 & idx <= length(TaskTimings));
        TaskTimings2(valid_idx) = 1;
    end

end

