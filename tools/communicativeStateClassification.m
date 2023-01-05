function [overlapWCH1, overlapWCH2, overlapBCH1, overlapBCH2, gapCH1, ...
  gapCH2, pauseCH1, pauseCH2, utteranceCH1, utteranceCH2, turnCH1, turnCH2] = ...
  communicativeStateClassification(actArr, binRes)
% SILENCEOVERLAPCLASSIFICATIONS
%    Classifies speech dominant portions 
%    
%    Inputs:
%    - actArr:   binary array indicating speech dominant (1)
%                or silence domaninant (0) windows.                         [2x* logical] 
%    - binRes:   time difference between two bins in actArr                 [1x1 double] 
%
%    Outputs:
%    - overlapW:          durations and idx of overlaps-within (overlaps
%                         within an interlocutor's speech stream)           [*x3 double]
%    - overlapB:          durations and idx of overlaps-between 
%                         (overlap between talkers during turn-taking)      [*x3 double]
%    - gap:               durations and idx of gaps (gap between talkers
%                         during turn-taking)                               [*x3 double]
%    - pause:             durations and idx of pauses (pauses in one  
%                         person's speech stream not interrupted  
%                         by the interlocutor)                              [*x3 double]
%    - utterance:         durations and idx of utterances (speech
%                         surrounded by 180 ms of silence excluding 
%                         overlaps-within                                   [*x3 double]

%% INPUT CHECK
if nargin < 2
  error('activity:argCheck','Too few input arguments');
elseif nargin > 2
  error('activity:argCheck','Too many input arguments');
end

%%
% Divide activity into two separate channels 
CH1 = actArr(1,:);
CH2 = actArr(2,:);

% Create leaveIdx struct holding scalar values for keeping track of when states are left
fields = {'idleCH1','idleCH2', 'activityCH1','activityCH2', ...
          'overlapCH2', 'overlapCH1'};
        
leaveIdx = struct(fields{1}, 0, fields{2}, 0, fields{3}, 0, ...
                  fields{4}, 0, fields{5}, 0, fields{6}, 0);
                
%% Initialization 
% Duration, start and stop indices arrays
overlapWCH1  = double.empty(0,3); % Overlap-within of CH1 in CH2's turn
overlapWCH2  = double.empty(0,3); % Overlap-within of CH2 in CH1's turn
overlapBCH1  = double.empty(0,3); % Overlap-between durations from CH2->CH1, so how fast CH1 responded
overlapBCH2  = double.empty(0,3); % Overlap-between durations from CH1->CH2, so how fast CH2 responded
gapCH1       = double.empty(0,3); % Gap durations from CH2->CH1, so how fast CH1 responded
gapCH2       = double.empty(0,3); % Gap durations from CH1->CH2, so how fast CH2 responded
pauseCH1     = double.empty(0,3); % Pause in CH1's turn
pauseCH2     = double.empty(0,3); % Pause in CH2's turn
utteranceCH1 = double.empty(0,3); % Interpausal units
utteranceCH2 = double.empty(0,3); % Interpausal units
turnCH1      = double.empty(0,3); % CH1's turn
turnCH2      = double.empty(0,3); % CH2's turn

% Initial state
state = s.idle;

%% Categorization
for i = 2:length(CH1)

  if or(xor(CH1(i), CH1(i-1)),xor(CH2(i), CH2(i-1))) % Detect any change

    % Combine binary to decimal states:
    leavePort = (sum([CH1(i), CH2(i)].*2.^(numel([CH1(i), CH2(i)])-1:-1:0)));
    % This will result in decoding of the following states:
    %  CH1  CH2
    %    0    0 => 0
    %    0    1 => 1
    %    1    0 => 2
    %    1    1 => 3  
    
    switch state
                
      case s.activityCH1        
        switch leavePort
          case 0 
            utteranceCH1(end+1,:) = getDurationAndIdx(i, [leaveIdx.idleCH1 leaveIdx.idleCH2 leaveIdx.activityCH2]);
            state = s.idleCH1;
          case 1 
            utteranceCH1(end+1,:) = getDurationAndIdx(i, [leaveIdx.idleCH1 leaveIdx.idleCH2 leaveIdx.activityCH2]);
            gapCH2(end+1,:) = getDurationAndIdx(i,i); % No-gap-no-overlap, i.e. gap = 0 ms
            turnCH1(end+1,:) = getDurationAndIdx(i, [leaveIdx.idleCH2 leaveIdx.activityCH2]);
            state = s.activityCH2;
          case 3
            state = s.overlapCH2;            
          otherwise
        end
        leaveIdx.activityCH1 = i; 
        
      case s.idleCH1
        switch leavePort
          case 1 
            gapCH2(end+1,:) = getDurationAndIdx(i, [leaveIdx.activityCH1 leaveIdx.overlapCH2]);
            turnCH1(end+1,:) = getDurationAndIdx([leaveIdx.activityCH1 leaveIdx.overlapCH2],...
              [leaveIdx.idleCH2 leaveIdx.activityCH2]);
            state = s.activityCH2;
          case 2 
            pauseCH1(end+1,:) = getDurationAndIdx(i, [leaveIdx.activityCH1 leaveIdx.overlapCH2]);
            state = s.activityCH1;
          case 3 
            pauseCH1(end+1,:) = getDurationAndIdx(i, [leaveIdx.activityCH1 leaveIdx.overlapCH2]);
            state = s.overlapCH2;
          otherwise 
        end
        leaveIdx.idleCH1 = i; 

      case s.overlapCH1
        switch leavePort
          case 0 
            utteranceCH2(end+1,:) = getDurationAndIdx(i, [leaveIdx.idleCH1 leaveIdx.idleCH2 leaveIdx.activityCH1]);
            overlapWCH1(end+1,:) = getDurationAndIdx(i, [leaveIdx.activityCH2 leaveIdx.idleCH2]);
            state = s.idleCH2;
          case 1 
            overlapWCH1(end+1,:) = getDurationAndIdx(i, [leaveIdx.activityCH2 leaveIdx.idleCH2]);
            state = s.activityCH2;
          case 2 
            utteranceCH2(end+1,:) = getDurationAndIdx(i, [leaveIdx.idleCH1 leaveIdx.idleCH2 leaveIdx.activityCH1]);
            overlapBCH1(end+1,:) = getDurationAndIdx(i, [leaveIdx.activityCH2 leaveIdx.idleCH2]);
            turnCH2(end+1,:)  = getDurationAndIdx(i, [leaveIdx.idleCH1 leaveIdx.activityCH1]);
            state = s.activityCH1;
          otherwise
        end
        leaveIdx.overlapCH1 = i;       

      case s.activityCH2   
        switch leavePort
          case 0 
            utteranceCH2(end+1,:) = getDurationAndIdx(i, [leaveIdx.idleCH1 leaveIdx.idleCH2 leaveIdx.activityCH1]);
            state = s.idleCH2;
          case 2 
            utteranceCH2(end+1,:) = getDurationAndIdx(i, [leaveIdx.idleCH1 leaveIdx.idleCH2 leaveIdx.activityCH1]);
            gapCH1(end+1,:) = getDurationAndIdx(i, i);  % No-gap-no-overlap, i.e. gap = 0 ms
            turnCH2(end+1,:) = getDurationAndIdx(i, [leaveIdx.idleCH1 leaveIdx.activityCH1]);
            state = s.activityCH1;
          case 3
            state = s.overlapCH1;
          otherwise
        end
        leaveIdx.activityCH2 = i; 
        
      case s.idleCH2
        switch leavePort
          case 1
            pauseCH2(end+1,:) = getDurationAndIdx(i, [leaveIdx.activityCH2 leaveIdx.overlapCH1]);
            state = s.activityCH2;
          case 2 
            gapCH1(end+1,:) = getDurationAndIdx(i, [leaveIdx.activityCH2 leaveIdx.overlapCH1]);
            turnCH2(end+1,:) = getDurationAndIdx([leaveIdx.activityCH2 leaveIdx.overlapCH1],...
              [leaveIdx.idleCH1 leaveIdx.activityCH1]);       
            state = s.activityCH1;
          case 3
            pauseCH2(end+1,:) = getDurationAndIdx(i, [leaveIdx.activityCH2 leaveIdx.overlapCH1]);
            state = s.overlapCH1;
          otherwise 
        end
        leaveIdx.idleCH2 = i; 
        
      case s.overlapCH2
        switch leavePort
          case 0
            utteranceCH1(end+1,:) = getDurationAndIdx(i, [leaveIdx.idleCH1 leaveIdx.idleCH2 leaveIdx.activityCH2]);
            overlapWCH2(end+1,:) = getDurationAndIdx(i, [leaveIdx.activityCH1 leaveIdx.idleCH1]);
            state = s.idleCH1;
          case 1 
            utteranceCH1(end+1,:) = getDurationAndIdx(i, [leaveIdx.idleCH1 leaveIdx.idleCH2 leaveIdx.activityCH2]);
            overlapBCH2(end+1,:) = getDurationAndIdx(i, [leaveIdx.activityCH1 leaveIdx.idleCH1]);
            turnCH1(end+1,:) = getDurationAndIdx(i, [leaveIdx.idleCH2 leaveIdx.activityCH2]);
            state = s.activityCH2;       
          case 2 
            overlapWCH2(end+1,:) = getDurationAndIdx(i, [leaveIdx.activityCH1 leaveIdx.idleCH1]);
            state = s.activityCH1;
          otherwise
        end
        leaveIdx.overlapCH2 = i;  
        
      otherwise % Initial state, loops until first activity
        switch leavePort
          case 1
            leaveIdx.idleCH2 = i; % Initial state seen as idleCH2 state
            state = s.activityCH2;
          case {2,3} % In case of competing first start, favor CH1 (choice)
            leaveIdx.idleCH1 = i; % Initial state seen as idleCH1 state
            state = s.activityCH1;
          otherwise
        end   
        
    end
  end
end

function returnArr = getDurationAndIdx(aEndIdx, aFromIdx)
  % aEndIdx  - Array containing possible "end" indices
  % aFromIdx - Array containing possible "from" indices
  
  % Determine final endIdx and fromIdx:
  endIdx  = max(aEndIdx);
  fromIdx = max(aFromIdx);
  % Calculate duration:
  duration = endIdx - fromIdx;
  
  returnArr = [duration fromIdx-1 endIdx-1]; % -1 cause MATLAB indices from 0
end

% Convert durations to seconds 
overlapWCH1(:,1)  = overlapWCH1(:,1)  * binRes;
overlapWCH2(:,1)  = overlapWCH2(:,1)  * binRes;
overlapBCH1(:,1)  = overlapBCH1(:,1)  * binRes;
overlapBCH2(:,1)  = overlapBCH2(:,1)  * binRes;
gapCH1(:,1)       = gapCH1(:,1)       * binRes;
gapCH2(:,1)       = gapCH2(:,1)       * binRes;
pauseCH1(:,1)     = pauseCH1(:,1)     * binRes;
utteranceCH1(:,1) = utteranceCH1(:,1) * binRes;
utteranceCH2(:,1) = utteranceCH2(:,1) * binRes;
turnCH1(:,1)      = turnCH1(:,1)      * binRes;
turnCH2(:,1)      = turnCH2(:,1)      * binRes;

end

