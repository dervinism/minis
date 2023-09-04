function validatedInput = validateInput(input)
% validatedInput = validateInput(input)
% Function validates input to minisHeadless function.
arguments
  input struct
end


%% Validate the main input part
p1 = inputParser;

addParameter(p1,'task','');

addParameter(p1,'loadTargetFileInput', '');
addParameter(p1,'loadNoiseFileInput', '');

addParameter(p1,'maxTimeToPeak', '10');
addParameter(p1,'baselineDuration', '2');
addParameter(p1,'peakIntegrationPeriod', '2.5');
addParameter(p1,'Amplobound', '0.02');
addParameter(p1,'Ampupbound', '10');
addParameter(p1,'smoothWindow', '1.5');
addParameter(p1,'RTint', '10-90%');
addParameter(p1,'downGoing', false);
addParameter(p1,'voltageClamp', false);

addParameter(p1,'RTbinSize', '0.5');

addParameter(p1,'startPulseTarget', '');
addParameter(p1,'endPulseTarget', '');
addParameter(p1,'startGlitchTarget', '');
addParameter(p1,'endGlitchTarget', '');
addParameter(p1,'startPulseNoise', '');
addParameter(p1,'endPulseNoise', '');
addParameter(p1,'startGlitchNoise', '');
addParameter(p1,'endGlitchNoise', '');
addParameter(p1,'pulseDuration', '0.5');

addParameter(p1,'loSimAmp', '0.01');
addParameter(p1,'L', '0.6');
addParameter(p1,'tau_m', '10');
addParameter(p1,'tau_PSPm', '20');

addParameter(p1,'distType', 'Quadrimodal Normal');
addParameter(p1,'distBaseline', 'Zero');
addParameter(p1,'SDupbound', '0.035');
addParameter(p1,'SDlobound', '0.025');
addParameter(p1,'maxErr', '2000');
addParameter(p1,'maxAmpErr', '1000');
addParameter(p1,'maxRTErr', '1000');
addParameter(p1,'maxAmpBottomErr', '1000');
addParameter(p1,'maxAmpMidErr', '1000');
addParameter(p1,'maxAmpTopErr', '400');
addParameter(p1,'maxDev', '150');
addParameter(p1,'maxDevAmp', '300');
addParameter(p1,'maxDevRT', '300');
addParameter(p1,'maxAmpBottomDev', '100');
addParameter(p1,'maxAmpMidDev', '100');
addParameter(p1,'maxAmpTopDev', '40');

addParameter(p1,'filtering', 'off');
addParameter(p1,'filtfs', '50, 150');

addParameter(p1,'options', '');

parse(p1,input)
input1 = rmfield(p1.Results, {'loadTargetFileInput','loadNoiseFileInput',...
  'downGoing','voltageClamp','options'});
input1 = [string(fieldnames(input1)), struct2cell(input1)].';
input1 = validate_p1("downGoing",input.downGoing,...
  "voltageClamp",input.voltageClamp,input1{:});


%% Validate target file
if ~ischar(p1.Results.loadTargetFileInput)
  p2 = inputParser;
  addParameter(p2,'dt', []);
  addParameter(p2,'lActualEpisodes', []);
  addParameter(p2,'sweep', []);
  addParameter(p2,'current', []);
  parse(p2,p1.Results.loadTargetFileInput)
  names = string(fieldnames(p2.Results));
  values = struct2cell(p2.Results);
  input1.loadTargetFileInput = validate_p2(names(1),values{1},...
    names(2),values{2},names(3),values{3},names(4),values{4});
else
  input1.loadTargetFileInput = p1.Results.loadTargetFileInput;
end


%% Validate noise file
if ~ischar(p1.Results.loadNoiseFileInput)
  p3 = inputParser;
  addParameter(p3,'dt', []);
  addParameter(p3,'lActualEpisodes', []);
  addParameter(p3,'sweep', []);
  addParameter(p3,'current', []);
  parse(p3,p1.Results.loadNoiseFileInput)
  names = string(fieldnames(p3.Results));
  values = struct2cell(p3.Results);
  input1.loadNoiseFileInput = validate_p2(names(1),values{1},...
    names(2),values{2},names(3),values{3},names(4),values{4});
else
  input1.loadNoiseFileInput = p1.Results.loadNoiseFileInput;
end


%% Validate options
p4 = inputParser;

addParameter(p4,'bounds', [.01,      .01,       0,   .25,      .01,    -1,   .01,      .01,  -10000,   .25,      .01,    -1,   .01,      .01,  -10000,   .25,      .01,    -1,     -10,     -10,     -10,     -10,     -10,     -10;
                             1,        1,   10000,    10,       10,     1,     1,        1,   10000,    10,       10,     1,     1,        1,   10000,    10,       10,     1,      10,      10,      10,      10,      10,      10]);
addParameter(p4,'nGenerations', 50);
addParameter(p4,'parallelCores', 'max');
addParameter(p4,'fullParallel', false);
addParameter(p4,'tauRange', false);
addParameter(p4,'cluster', false);
addParameter(p4,'clusterProfile', 'local');
addParameter(p4,'cliff', false);
addParameter(p4,'figureDisplay', true);
addParameter(p4,'SDlobound', false);
addParameter(p4,'SDupbound', false);
addParameter(p4,'estimateTauLo', false);
addParameter(p4,'estimateTauHi', false);
addParameter(p4,'optimisationData', '');
addParameter(p4,'resumeOptimisation', '');

parse(p4,input.options)
names = string(fieldnames(p4.Results));
values = struct2cell(p4.Results);
input2 = validate_p4(names(1),values{1},names(2),values{2},...
  names(3),values{3},names(4),values{4},names(5),values{5},...
  names(6),values{6},names(7),values{7},names(8),values{8},...
  names(9),values{9},names(10),values{10},names(11),values{11},...
  names(12),values{12},names(13),values{13},names(14),values{14},...
  names(15),values{15});


%% Concatenate validated input
validatedInput.task = input1.task;
validatedInput.loadTargetFileInput = input1.loadTargetFileInput;
validatedInput.loadNoiseFileInput = input1.loadNoiseFileInput;
validatedInput.maxTimeToPeakEdit = input1.maxTimeToPeak;
validatedInput.baselineDurationEdit = input1.baselineDuration;
validatedInput.peakIntegrationPeriodEdit = input1.peakIntegrationPeriod;
validatedInput.AmploboundEdit = input1.Amplobound;
validatedInput.AmpupboundEdit = input1.Ampupbound;
validatedInput.smoothWindowEdit = input1.smoothWindow;
validatedInput.RTintEdit = input1.RTint;
validatedInput.downGoingCheckbox = input1.downGoing;
validatedInput.voltageClampCheckbox = input1.voltageClamp;
validatedInput.RTbinSizeEdit = input1.RTbinSize;
validatedInput.startPulseTargetEdit = input1.startPulseTarget;
validatedInput.endPulseTargetEdit = input1.endPulseTarget;
validatedInput.startGlitchTargetEdit = input1.startGlitchTarget;
validatedInput.endGlitchTargetEdit = input1.endGlitchTarget;
validatedInput.startPulseNoiseEdit = input1.startPulseNoise;
validatedInput.endPulseNoiseEdit = input1.endPulseNoise;
validatedInput.startGlitchNoiseEdit = input1.startGlitchNoise;
validatedInput.endGlitchNoiseEdit = input1.endGlitchNoise;
validatedInput.pulseDurationEdit = input1.pulseDuration;
validatedInput.loSimAmpEdit = input1.loSimAmp;
validatedInput.LEdit = input1.L;
validatedInput.tau_mEdit = input1.tau_m;
validatedInput.tau_PSPmEdit = input1.tau_PSPm;
validatedInput.distTypeEdit = input1.distType;
validatedInput.distBaselineEdit = input1.distBaseline;
validatedInput.SDupbound = input1.SDupbound;
validatedInput.SDlobound = input1.SDlobound;
validatedInput.maxErrEdit = input1.maxErr;
validatedInput.maxAmpErrEdit = input1.maxAmpErr;
validatedInput.maxRTErrEdit = input1.maxRTErr;
validatedInput.maxAmpBottomErrEdit = input1.maxAmpBottomErr;
validatedInput.maxAmpMidErrEdit = input1.maxAmpMidErr;
validatedInput.maxAmpTopErrEdit = input1.maxAmpTopErr;
validatedInput.maxDevEdit = input1.maxDev;
validatedInput.maxDevAmpEdit = input1.maxDevAmp;
validatedInput.maxDevRTEdit = input1.maxDevRT;
validatedInput.maxAmpBottomDevEdit = input1.maxAmpBottomDev;
validatedInput.maxAmpMidDevEdit = input1.maxAmpMidDev;
validatedInput.maxAmpTopDevEdit = input1.maxAmpTopDev;
validatedInput.filtering = input1.filtering;
validatedInput.filtfs = {input1.filtfs};
validatedInput.options = input2;



%% Local functions
function input = validate_p1(input)
arguments
  input.task char {mustBeMember(input.task,["preprocess","detection",...
    "detectionHeadless","detectCompare","errorBounds","autoDistributionFit",...
    "autoDistributionFitHeadless","simulation","simulationHeadless"])}

  input.maxTimeToPeak char = '10'
  input.baselineDuration char = '2'
  input.peakIntegrationPeriod char = '2.5'
  input.Amplobound char = '0.02'
  input.Ampupbound char = '10'
  input.smoothWindow char = '1.5'
  input.RTint char {mustBeMember(input.RTint,["10-90%","20-80%"])} = '10-90%'
  input.downGoing logical = false
  input.voltageClamp logical = false

  input.RTbinSize char {mustBeMember(input.RTbinSize,["0.25","0.5"])} = '0.5'

  input.startPulseTarget char = ''
  input.endPulseTarget char = ''
  input.startGlitchTarget char = ''
  input.endGlitchTarget char = ''
  input.startPulseNoise char = ''
  input.endPulseNoise char = ''
  input.startGlitchNoise char = ''
  input.endGlitchNoise char = ''
  input.pulseDuration char = '0.5'

  input.loSimAmp char = '0.01'
  input.L char = '0.6'
  input.tau_m char = '10'
  input.tau_PSPm char = '20'

  input.distType char {mustBeMember(input.distType,["Normal",...
    "Bimodal Normal","Trimodal Normal","Quadrimodal Normal","Log Normal",...
    "Bimodal Log Normal","Trimodal Log Normal","Gaussian",...
    "Bimodal Gaussian","Trimodal Gaussian","Quadrimodal Gaussian",...
    "Skew-normal","Bimodal Skew-normal","Trimodal Skew-normal"])} = 'Quadrimodal Normal'
  input.distBaseline char {mustBeMember(input.distBaseline,["Zero",...
    "Subtracted"])} = 'Zero'
  input.SDupbound char = '0.035'
  input.SDlobound char = '0.025'
  input.maxErr char = '2000'
  input.maxAmpErr char = '1000'
  input.maxRTErr char = '1000'
  input.maxAmpBottomErr char = '1000'
  input.maxAmpMidErr char = '1000'
  input.maxAmpTopErr char = '400'
  input.maxDev char = '150'
  input.maxDevAmp char = '300'
  input.maxDevRT char = '300'
  input.maxAmpBottomDev char = '100'
  input.maxAmpMidDev char = '100'
  input.maxAmpTopDev char = '40'
  
  input.filtering char = 'off'
  input.filtfs char = '50, 150'
end

function input = validate_p2(input)
arguments
  input.dt (1,1) double
  input.lActualEpisodes (1,1) double
  input.sweep (1,:) double
  input.current (1,:) double
end

function input = validate_p4(input)
arguments
  input.bounds (2,:) {mustBeNumeric,mustBeReal} = [.01,      .01,       0,   .25,      .01,    -1,   .01,      .01,  -10000,   .25,      .01,    -1,   .01,      .01,  -10000,   .25,      .01,    -1,     -10,     -10,     -10,     -10,     -10,     -10;
                                                     1,        1,   10000,    10,       10,     1,     1,        1,   10000,    10,       10,     1,     1,        1,   10000,    10,       10,     1,      10,      10,      10,      10,      10,      10];
  input.nGenerations (1,1) {mustBeNumeric,mustBeReal,mustBePositive} = 50
  input.parallelCores char = 'max'
  input.fullParallel logical = false
  input.tauRange logical = false
  input.cluster logical = false
  input.clusterProfile char = 'local'
  input.cliff logical = false
  input.figureDisplay logical = true
  input.SDlobound logical = false
  input.SDupbound logical = false
  input.estimateTauLo logical = false
  input.estimateTauHi logical = false
  input.optimisationData char = ''
  input.resumeOptimisation char = ''
end