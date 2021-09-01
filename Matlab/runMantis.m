function varargout = runMantis(Scenario, client, varargin)
%UNTITLED5 runMantis run the simulation for the options specified in the
%scenario.
%   This function triggers a sequence of events
%   1 prints the scenario options to a file.
%   2 Executes the test clinet program which reads the input file 
%     converts it to string message and sends it to the mantisServer
%   3 MantisServer runs the simulation and sends the results back to the 
%     client program
%   4 The client program receives the message and prints the results to a
%     file
%   5 Last this script reads the file and returns the results to the matlab
%     workspace
%
%   Scenario:   is a structure with the scenario options. You can get one from
%               the MantisInputs
%   client:     is the paths including the executable name of the test client program      

if ~isempty(varargin)
   if strcmp('quit',  varargin{1})
       system([client ' quit']);
       return
   end
end

delete(Scenario.outfile)
writeMantisInput(Scenario);
system([client ' ' Scenario.infile ' ' Scenario.outfile]);
[btc tf] = readMantisOutput(Scenario.outfile);
varargout{1} = btc;
varargout{2} = tf;


end

