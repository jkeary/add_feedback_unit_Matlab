% Name: James Keary 
% Student ID: N12432851 
% Net ID: jpk349 
%
% Class: MPATE-GE.2599 - Fundamentals of Digital Signal Theory
% Semester: Fall 2011
% Instructor:  Dr. Agnieszk Roginska
%
% Assignment: DST Final Project.
% Due: 12/19/2011 @ 11:59 pm
%
% ---------- FUNCTION DESCRIPTION ----------
%
% addfeedback ( INfilename, partial, attack, decay, sustain,
%               release, OUTfilename, fundFreq )
%
% The function will add digital feedback to an input file based on
% parameters input by the user.  The function's output is a file of the
% input file with feedback at specified user-input parameters.
%
% Inputs
%
%   INfilename:     Name of .wav file containing audio of a plucked guitar
%                   string
%
%   partial:        the partial number the user wants to feedback:
%                       1 = 1st partial: fundamental frequency
%                       2 = 2nd partial: 1st overtone, octave above, 
%                           fundamental x 2 
%                       3 = 3rd partial: 2nd overtone, fifth above,
%                           fundamental x 3
%
%   attack:         attack of the feedback frequency envelope in seconds
%
%   decay:          decay of the feedback frequency envelope in seconds
%
%   sustain:        sustain of the feedback frequency envelope in seconds
%
%   release:        release of the feedback frequency envelope in seconds
%
%   OUTfilename:    Name of .wav file to which the resulting (feedback)
%                   signal will be written   
%
%   %%An Optional input%%
%
%   fundFreq:       User has the option to input the fundamental frequency
%                   of the .wav file.  Here is a reference table of the
%                   fundamental frequencies (Hz) of the example .wav files
%                   provided with the function:
%                       'lowE string.wav'   : fundFreq = 82.41;
%                       'A string.wav'      : fundFreq = 110.00;
%                       'D string.wav'      : fundFreq = 146.83;
%                       'G string.wav'      : fundFreq = 196.00;
%                       'B string.wav'      : fundFreq = 246.94;
%                       'HiE string.wav'    : fundFreq = 329.63;
%                   IF THE USER DECIDES TO INPUT THEIR OWN .WAV FILE, THEN
%                   THEY ALSO SHOULD ENTER THE fundFreq INPUT.
%
% Output
%
%   A file if the INFile with synthetic feedback added based on user input

function addfeedback ( INfilename, partial, attack, decay, sustain, release, OUTfilename, fundFreq )

%-------- ERROR CHECKING OF INPUTS AND ASSIGNING VARIABLES ---------

% char checks: make sure your filenames are strings.

if ischar(INfilename) == 0,
    error('INfilename input needs to be a string');
end

% number of arguments check.

if nargin < 7,
    error('Not enough inputs, need at least 7')
end

if nargin > 8,
    error('Too many inputs, need at most 8')
end

% partial check

if partial > 3 || partial < 0,
    error('partial must be 1, 2, or 3')
end

% amplitude envelope of the feedback check

if isnumeric(attack) == 0 || isnumeric(decay) == 0 || isnumeric(sustain) == 0 || isnumeric(release) == 0,
    error('all of your feedback envelope parameters must be numbers')
end


% ---------- FUNCTION COMPUTATIONS ----------

% The function reads the .WAV file, converts to mono and normalizes to
% prevent clipping
[y, fs] = wavread(INfilename);
yMono = mean(y,2);
yNormal = yMono / 1.001 * max(abs(yMono));

% Function assigns frequencies to input file
if nargin == 7
    switch INfilename
        case 'lowE string.wav'
            fundFreq = 82.41;
        case 'A string.wav'
            fundFreq = 110.00;
        case 'D string.wav'
            fundFreq = 146.83;
        case 'G string.wav'
            fundFreq = 196.00;
        case 'B string.wav'
            fundFreq = 246.94;
        case 'HiE string.wav'
            fundFreq = 329.63;
        otherwise
            error('You have either incorrectly entered one of the example .wav files I provided with the function, or you are attempting to input your own .wav file.  If you are entering your own, please provide the fundFreq, input argument 8.')
    end
end

if nargin == 8
    switch INfilename
        case 'lowE string.wav'
            fundFreq = 82.41;
        case 'A string.wav'
            fundFreq = 110.00;
        case 'D string.wav'
            fundFreq = 146.83;
        case 'G string.wav'
            fundFreq = 196.00;
        case 'B string.wav'
            fundFreq = 246.94;
        case 'HiE string.wav'
            fundFreq = 329.63;
        otherwise
    end
end

% Function assigns feedback frequency based on user input of partial
switch partial
    case 1
        feedFreq = fundFreq;
    case 2
        feedFreq = 2*fundFreq;
    case 3
        feedFreq = 3*fundFreq;
end

% Function determines the number of samples
numSamps = ceil(fs / feedFreq)*2;

% Function generates filter coefficients
[b,a] = butter(1, .9, 'low');

% Function takes a section of samples from the .wav file, something
% after the initial transients of the plucked string have died down.
section = yNormal(39500:39500+numSamps, :);
% plot(section)

% Function creates ADSR envelope from user inputs
atk = linspace(0,1, attack*fs);
dk = linspace(1,.5, decay*fs);
sus = linspace(.5, .6, sustain*fs);
rel = linspace(.6, 0, release*fs);
ADSR = [atk,dk,sus,rel]';

% Function zero-pads
if length(ADSR) > length(yNormal),
    zeropad = length(ADSR)-length(yNormal);
    yNormal = [yNormal; zeros(zeropad,1)];
end

if length(ADSR) < length(yNormal),
    zeropad = length(yNormal)-length(ADSR);
    ADSR = [ADSR; zeros(zeropad,1)]; 
end

% Function creates hops and other variables for for loop.  
hop = floor(numSamps/2);
loops = floor(length(ADSR)/hop-1);
triwin = hann(numSamps);
outputMatrix = zeros(length(ADSR), 1);
start = 1; finish = numSamps;

% Function overlap adds section looping ontop of itself while filtering
% frequencies and windowing the section to prevent clipping.  The looping
% of the filter creates changes in the attenuation of the phase of the
% section, so the loop has phase movement over time.
for k = 1:loops
    section = (filter(b, a, section));
    temp = section(1:numSamps)* 1000000  .* triwin ; 
    outputMatrix(start:finish, :) = outputMatrix(start:finish, :) + temp; 
    start = start + hop;
    finish = finish + hop;
end

% function puts together ADSR and outputMatrix
outputMatrix = (outputMatrix ./ max(abs(outputMatrix))) .* ADSR;
outputMatrix = outputMatrix + yNormal; 
outputMatrix = outputMatrix ./ max(abs(outputMatrix));
% plot(outputMatrix);

% ---------- FUNCTION OUTPUT ----------

% Function creates file
wavwrite (outputMatrix, fs, OUTfilename);

end

