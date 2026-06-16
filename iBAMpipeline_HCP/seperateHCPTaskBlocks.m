function [neutTaskTimings,fearTaskTimings]=seperateHCPTaskBlocks(TaskTimings)
%   DESCRIPTION:
%       This function serves as a preprocessing utility for the iBAM framework 
%       to isolate experimental conditions within the Human Connectome Project 
%       (HCP) Emotion Processing Task. It splits the unified task timing vector 
%       into two distinct paradigm structures: 'Neutral' (shape matching) and 
%       'Fear' (face matching) blocks.
%
%       CRITICAL TEMPORAL GUARDRAILING:
%       Given the highly compact design of the HCP task blocks and the rapid 
%       alternation between conditions, the blood-oxygen-level-dependent (BOLD) 
%       hemodynamic response or sliding-window features are highly susceptible 
%       to cross-condition contamination. 
%
%       To mitigate temporal spillover effects, this function automatically 
%       truncates the boundary frames. By default, it discards the first 5 frames 
%       and the last 5 frames of each isolated task block, ensuring maximum statistical purity for downstream mapping.

neutTaskTimings=zeros(size(TaskTimings));
neutTaskTimings(TaskTimings==-1)=1;
for i=5:length(TaskTimings)-5
    if TaskTimings(i)==0&& TaskTimings(i+1)==-1
       neutTaskTimings(i+1:i+5)=0;
    end

end
for i=5:length(TaskTimings)-5
    if TaskTimings(i)==0&& TaskTimings(i-1)==-1
       neutTaskTimings(i-5:i-1)=0;
    end

end


fearTaskTimings=zeros(size(TaskTimings));
fearTaskTimings(TaskTimings==1)=1;
for i=5:length(TaskTimings)-5
    if TaskTimings(i)==0&& TaskTimings(i+1)==1
       fearTaskTimings(i+1:i+5)=0;
    end

end
for i=5:length(TaskTimings)-5
    if TaskTimings(i)==0&& TaskTimings(i-1)==1
       fearTaskTimings(i-5:i-1)=0;
    end

end