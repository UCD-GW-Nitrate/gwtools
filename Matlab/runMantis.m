function varargout = runMantis(Scenario)
% [btc, tf] = runMantis(Scenario) run the simulation for the options specified in the scenario.
%   This function triggers a sequence of events
%   1 prints the scenario options to a file.
%   2 Executes the test client program which reads the input file, 
%     converts it to string message, and sends it to the mantisServer
%   3 MantisServer runs the simulation and sends the results back to the 
%     client program
%   4 The client program receives the message and prints the results to a
%     file
%   5 Last this script reads the file and returns the results to the matlab
%     workspace
%
%   Scenario:   is a structure with the scenario options. You can get one
%               by running MantisInputs()
 
if isempty(Scenario.client)
    warning('Scenario.client is empty');
    for ii = 1:nargout
        varargout{ii} = [];
    end
    return
end

% if ~isempty(varargin)
%    if strcmp('quit',  varargin{1})
%        system([Scenario.client ' quit > $null']);
%        return
%    end
% end

delete(Scenario.outfile);
writeMantisInput(Scenario);
system([Scenario.client ' ' Scenario.infile ' ' Scenario.outfile '> $null']);
[btc, tf] = readMantisOutput(Scenario.outfile);
if nargout > 0
    varargout{1} = btc;
end
if nargout > 1
    varargout{2} = tf;
end


end

