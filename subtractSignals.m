function [ data ] = subtractSignals( filename, linAndSineFitDomain )
% subtractSignals imports an XLS or XLSX file and subtracts both a linear
% background and a sum of Sines fringe signal
%   Detailed explanation goes here
disp('Importing data');
sw = warning; % save warning state
try
    warning('error', 'stats:dataset:genvalidnames:ModifiedVarnames'); % convert warning to error
    fromFile = dataset('XLSFile', filename); % import XLS as dataset
catch
    warning('off','stats:dataset:genvalidnames:ModifiedVarnames'); % turn off warning
    fromFile = dataset('XLSFile', filename); % import XLS as dataset
    % format Variable Names for printing
    varNames = fromFile.Properties.VarNames{1};
    for i=2:length(fromFile.Properties.VarNames)
        varNames = [varNames ', ' fromFile.Properties.VarNames{i}];
    end
    warning('Variable names were modified to make them valid MATLAB identifiers. Variable Names are now %s', varNames);
end   
warning(sw); % restore warning state
disp('Cropping data for background and fringe estimation');
dataToFit = cropData(fromFile, linAndSineFitDomain); % crop data for fitting background/line, and fringes/Sine
disp('Estimating background');
lineCoeffs = fitLineToData(dataToFit); % fit linear model func to cropped data
disp('Estimating fringes');
dataToFitMinusLine = subtractLine(dataToFit, lineCoeffs); % subtract linear func from cropped data
sineCoeffs = fitSineToData(dataToFitMinusLine); % fit sum of Sines model func to cropped data minus line
disp('Subtracting background and fringe signals');
%rotate all curves by first curve's line fit 
for i=2:length(lineCoeffs.Properties.VarNames)
    lineCoeffs(:,i)=lineCoeffs(:,1);
end
dataMinusLine = subtractLine(fromFile, lineCoeffs); % subtract linear func from data
dataMinusSine = subtractSines(fromFile, sineCoeffs); % subtract sum of Sines func from data
dataMinusSineMinusLine = subtractLine(dataMinusSine, lineCoeffs); % subtract linear func from data less sum of Sines func
disp('Consolidating data');
data = horzcat(fromFile,dataMinusLine(:,2:N), dataMinusSine(:,2:N), dataMinusSineMinusLine(:,2:N)); % concatenate datasets without repeating x/WaveNumber

% generate plots tiled 2x2
disp('Generating plot');
N = length(fromFile.Properties.VarNames); % get number of columns now because MATLAB has stupid scoping rules
% top-left plot is original data
subplot(2,2,1), plot(fromFile(:,1),fromFile(:,2:N)), 
xlabel(fromFile.Properties.VarNames{1}), ylabel('intensity'), 
title('Original Data'), legend(fromFile.Properties.VarNames(2:N)); 

% top-right plot is original data minus background(line)
subplot(2,2,2), plot(dataMinusLine(:,1),dataMinusLine(:,2:N)), 
xlabel(fromFile.Properties.VarNames{1}), ylabel('intensity'), 
title('Original Minus BackGround(Line)');

% bottom-left plot is original data minus fringes (sum of Sines)
subplot(2,2,3), plot(dataMinusSine(:,1),dataMinusSine(:,2:N)), 
xlabel(fromFile.Properties.VarNames{1}), ylabel('intensity'), 
title('Original Minus Fringes(Sine)');

% bottom-right plot is original data minus both background(line) and
% fringes(sum of Sines)
subplot(2,2,4), plot(dataMinusSineMinusLine(:,1),dataMinusSineMinusLine(:,2:N)), 
xlabel(fromFile.Properties.VarNames{1}), ylabel('intensity'), 
title('Original Minus BackGround(Line) and Fringes(Sine)');

% Write to file
writeToFile = input('Write data to CSV file? y/n [n]: ','s'); % prompt user to save data
if isempty(writeToFile)
    writeToFile = 'n'; % default answer is no
end
if writeToFile=='y'
    filename = input('Enter a filename [signal_subtraction.csv]: ', 's'); % prompt user for filename 
    if isempty(filename)
        filename = 'signal_subtraction.csv'; % default filename is signal_subtraction.csv
    end
    disp('Writing to file');
    export(data,'File',filename,'delimiter',','); % write data to CSV file
end
disp('Goodbye');

% private functions
function [ croppedData ] = cropData( data, domainRange )
%cropData crops data to include only x and f_1(x),...f_N(x) s.t.
%domainRange(1) <= x <= domainRange(2)
%   Detailed explanation goes here
i = find(data.(data.Properties.VarNames{1})>=domainRange(1), 1 ); % get smallest index of x such that x >= domainRange(1)
j = find(data.(data.Properties.VarNames{1})<=domainRange(2), 1, 'last'); % get largest index of x such that x <= domainRange(2)
croppedData = data(i:j,:); % crop data, keeping data(k,:), i<=k<=j
end
function [ coeffs ] = fitLineToData( data )
%fitLineToData fits a linear model function to x,yi pairs in data
%   Detailed explanation goes here
N = length(data.Properties.VarNames); % get number of columns
x = data.(data.Properties.VarNames{1}); % get x values
c = zeros(2,N-1); % preallocate coefficient matrix
% fit line to y1,y2,...,yN-1
for i=2:N
    y = data.(data.Properties.VarNames{i}); % get y values
    c(:,i-1) = polyfit(x, y, 1); % call polyfit on x,y and store [m;x] in column i-1 of c
end
coeffs = mat2dataset(c,'VarNames',data.Properties.VarNames(2:N));
end
function [ coeffs ] = fitSineToData( data )
% fitSineToData fits a sum of Sines model function to x,yi pairs in data
%   Details
N = length(data.Properties.VarNames); % get number of columns
x = data.(data.Properties.VarNames{1}); % get x values
c = zeros(6,N-1); % preallocate coefficient matrix
% fit line to y1,y2,...,yN-1
for i=2:N
    y = data.(data.Properties.VarNames{i}); % get y values
    f = fit(x, y, 'sin2'); % fit sum of Sines model of size 2 to data
    c(:,i-1) = coeffvalues(f)'; % store a1, ..., c2 in column i-1 of c
end
coeffs = mat2dataset(c,'VarNames',data.Properties.VarNames(2:N));
end
function [ dataMinusLine ] = subtractLine( data, coeffs)
%subtractLine subtracts a line from each signal in data
%   Detailed explanation goes here
M = length(data.(data.Properties.VarNames{1})); % get number of rows
x = data.(data.Properties.VarNames{1}); % store wavenumber as x
lines = [x, ones(M,1)]*double(coeffs); % calculate lines for each set of coeffs
newVarNames = strcat(data.Properties.VarNames, '_rotated'); % concatenate all VarNames with ' - line'
dataMinusLine = mat2dataset(double(data) - [zeros(M,1),lines],'VarNames', newVarNames); % calculate signals - lines and store in a dataset 
end
function [ dataMinusSines ] = subtractSines( data, coeffs )
%subtractSines subtracts a sum of Sines from signals in data given a
%coefficient matrix
%   Detailed explanation goes here
M = length(data.(data.Properties.VarNames{1})); % get number of rows
x = data.(data.Properties.VarNames{1}); % store wavenumber as x
N = length(data.Properties.VarNames)-1; % N given y_1,y_2,...,y_N
sumOfSines = zeros(M,N); % preallocate sumOfSines 
sin2 = fittype('sin2'); % get fittype object
% calculate sum of Sines
for i = 1:N
    c = coeffs.(coeffs.Properties.VarNames{i}); % get coeffs for y_i
    sumOfSines(:,i) = feval(sin2,c(1),c(2),c(3),c(4),c(5),c(6),x); % evaluate for y_i
end
newVarNames = strcat(data.Properties.VarNames, '_fringesRemoved'); % concatenate all VarNames with '_fringesRemoved'
dataMinusSines = mat2dataset(double(data) - [zeros(M,1),sumOfSines],'VarNames', newVarNames); % calculate signals - Sines and store in a dataset 
end
end

