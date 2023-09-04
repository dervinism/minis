function varargout = minisInit(varargin)
% minisInit
% [status, output] = minisInit(input);
%
% Function initialises minis software. If supplied without any input,
% initialises minis GUI and does not produce any output variables. If
% supplied with an input structure variable (type help testMatlab for an
% example of how to set up the input structure variable), runs minis
% software without the GUI. Output variable 'status' reports the state of
% the algorithm execution with 0 for success and -1 otherwise. 'ouput' is a
% structure variable with fields explained in the minisHeadless function.
% For the definition of 'output' type help minisHeadless.

if nargin == 0
  minis;
elseif nargin == 1
  [varargout{1}, varargout{2}] = minisHeadless(varargin{1});
end