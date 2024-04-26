function id = ICHNOSExitReason2Num(ER)
%id = ICHNOSExitReason2Num(ER) Converts the exit reason to code
% 1 -> EXIT_TOP
% 2 -> EXIT_SIDE
% 3 -> EXIT_BOTTOM
% 4 -> MAX_INNER_ITER
% 5 -> STUCK
% 6 -> MAX_AGE
% 7 -> ATTRACT
% 8 -> INIT_OUT
% 9 -> FIRST_POINT_GHOST
% 10 -> FAR_AWAY
% 11 -> NOT_IN_SUBDOMAIN
% 12 -> CHANGE_PROCESSOR
% 13 -> NO_EXIT
% 14 -> NO_REASON


if strcmp(ER,'EXIT_TOP')
    id = 1;
elseif strcmp(ER,'EXIT_SIDE')
    id = 2;
elseif strcmp(ER,'EXIT_BOTTOM')
    id = 3;
elseif strcmp(ER,'MAX_INNER_ITER')
    id = 4;
elseif strcmp(ER,'STUCK')
    id = 5;
elseif strcmp(ER,'MAX_AGE')
    id = 6;
elseif strcmp(ER,'ATTRACT')
    id = 7;
elseif strcmp(ER,'INIT_OUT')
    id = 8;
elseif strcmp(ER,'FIRST_POINT_GHOST')
    id = 9;
elseif strcmp(ER,'FAR_AWAY')
    id = 10;
elseif strcmp(ER,'NOT_IN_SUBDOMAIN')
    id = 11;
elseif strcmp(ER,'CHANGE_PROCESSOR')
    id = 12;
elseif strcmp(ER,'NO_EXIT')
    id = 13;
elseif strcmp(ER,'NO_REASON')
    id = 14;
end