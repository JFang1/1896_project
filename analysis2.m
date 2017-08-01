%%--------------------------------------------------------%%
%
%Name: analysis2.m
%
%Description: This script is intended to be used with the 
%             wearable step tracking device developed as a 
%             senior design project during the Summer 2017
%             semester at the University of Pittsburgh. It
%             It uses the data acquired from the device
%             to estimate stride and step distances
%             during a straight line walk test.
%
%Authors:
%   Nathan Anuskiewicz
%   Jennifer Fang
%   Michael Mortensen
%
%   Last Modified:      Modified By:            Comments:
%   7/31/2017           Nathan Anuskiewicz      Initial Version
%
%%--------------------------------------------------------%%


clc

% request overall test length from user
len = input('Enter length of test (in meters): ');

% importing data for right foot
rightFile = 'RIGHT.txt';
[rAccel,rDelimeterOut] = importdata(rightFile);

% importing data for left foot
leftFile = 'LEFT.txt';
[lAccel,lDelimeterOut] = importdata(leftFile);

% initializing arrays that will hold velocity and displacement data
%  for the right foot
rV = zeros(size(rAccel,1),3); % 3 columns, for x/y/z
rD = zeros(size(rAccel,1),3);
%  and for the left foot
lV = zeros(size(lAccel,1),3);
lD = zeros(size(lAccel,1),3);
% initialize array for time intervals
rT = zeros(size(rAccel,1),1);
lT = zeros(size(lAccel,1),1);

% plot original accel
% figure(1);
% plot(lAccel(:,4),lAccel(:,1));
% title('original left accel');

% define new array of X Y Z accel
ravg_x = rAccel(:,1);
ravg_y = rAccel(:,2);
ravg_z = rAccel(:,3);

lavg_x = lAccel(:,1);
lavg_y = lAccel(:,2);
lavg_z = lAccel(:,3);

% creating a smoothed acceleration curve using rolling average
movAvg = 14;
coeff = ones(1,movAvg)/movAvg;

% smooth right foot curve
avg_x = filter(coeff,1,ravg_x);
avg_y = filter(coeff,1,ravg_y);
avg_z = filter(coeff,1,ravg_z);
smoothAccelR = [avg_x,avg_y,avg_z, rAccel(:,4)];

% smooth left foot curve
avg_x = filter(coeff,1,lavg_x);
avg_y = filter(coeff,1,lavg_y);
avg_z = filter(coeff,1,lavg_z);
smoothAccelL = [avg_x,avg_y,avg_z, lAccel(:,4)];

% calculating velocity data for the right foot
rAccelMag = abs(smoothAccelR);
rHeelStrikes = rAccelMag < .8; % TODO: might have to change this
rHeelLift = rAccelMag < .6;
% previous state (0 = stationary, 1 = moving)
prevState = 0;
smoothAccelR(:,1) = smoothAccelR(:,1)*-1;
% disp('debugging lines 84-101');
for w = 3:length(rV)-2
    % get time displacement for right foot
    rT(w)= smoothAccelR(w,4) - smoothAccelR(w-1,4);
    if(rHeelLift(w-2) == 0 && rHeelLift(w+2)==0 && rHeelLift(w-1) == 0 && rHeelLift(w)==0 && rHeelLift(w+1)==0 || prevState==1)
        rV(w,:) = rV(w-1,:) + smoothAccelR(w,1:3) * rT(w)*(.001);
        % zero Out any backwards velocities
        if( rV(w,1) < 0)
            rV(w,:)=[0,0,0];
        end
        prevState = 1;
    end
    if(rHeelStrikes(w-2) == 1 && rHeelStrikes(w+2) == 1 && rHeelStrikes(w-1) == 1 && rHeelStrikes(w) == 1 && rHeelStrikes(w+1) == 1 || prevState==0)
        rV(w,:) = [0 0 0]; % force zero velocity when foot stationary
        prevState = 0;
    end
end

% calculating velocity data for the left foot
lAccelMag = abs(smoothAccelL);
lHeelStrikes = lAccelMag < .8;
lHeelLift = lAccelMag < .6;
prevState=0;
for w = 3:length(lV)-2
    % get time displacement for right foot
    lT(w)= smoothAccelL(w,4) - smoothAccelL(w-1,4);
    if(lHeelLift(w-2) == 0 && lHeelLift(w+2) == 0 && lHeelLift(w-1) == 0 && lHeelLift(w) == 0 && lHeelLift(w+1) == 0 || prevState==1)
        lV(w,:) = lV(w-1,:) + smoothAccelL(w,1:3) * lT(w)*(.001);
        % zero Out any backwards velocities
         if( lV(w,1) < 0)
             lV(w,:)=[0,0,0];
         end
        prevState=1;
    end
    if(lHeelStrikes(w-2) == 1 && lHeelStrikes(w+2) == 1 && lHeelStrikes(w-1) == 1 && lHeelStrikes(w) == 1 && lHeelStrikes(w+1) == 1 || prevState==0)
        lV(w,:) = [0 0 0]; % force zero velocity when foot stationary
        prevState=0;
    end
end

% smooth velocity
movAvgv = 8;
coeffv = ones(1,movAvgv)/movAvgv;
rV = filter(coeffv,1,rV);
lV = filter(coeffv,1,lV);

% calculating displacement for the right foot
for ri = 2:size(smoothAccelR,1)
    for ry = 1:3
        rD(ri,ry) = rD(ri-1,ry) + rV(ri,ry) * rT(ri)*(.001);
    end
end

% calculating displacement for the left foot
for li = 2:size(smoothAccelL,1)
    for ly = 1:3
        lD(li,ly) = lD(li-1,ly) + lV(li,ly) * lT(w)*(.001);
    end
end

% smooth displacement
movAvgd = 6;
coeffd = ones(1,movAvgd)/movAvgd;
rD = filter(coeffd,1,rD);
lD = filter(coeffd,1,lD);

% determining individual stride lengths . . .

% array marker for right and left heel strike data in x+z dimensions
%  other calculation results.
rStrideD = zeros(size(smoothAccelR,1), 2);
lStrideD = zeros(size(smoothAccelL,1), 2);
% iterators
rj = 1;
lj = 1;

% finding the displacements at heel strike - right foot
for ri = 1:size(rStrideD,1) % for every sample
    if rV(ri,1) == 0 && rV(ri,2) == 0 && rV(ri,3) == 0 % if heelstrike
        rStrideD(rj,1) = rD(ri,1); % x dimension of displacement
        rj = rj + 1; % increment row in rStrideD
    end
end
% finding the displacements at heel strike - left foot
for li = 1:size(lStrideD,1) % for every sample
    if lV(li,1) == 0 && lV(li,2) == 0 && lV(li,3) == 0 % if heelstrike
        lStrideD(lj,1) = lD(li,1); % x dimension of displacement
        lj = lj + 1; % increment row in lStrideD
    end
end

% getting rid of repetition in the stride data
rStrideDExtracted = zeros(floor(size(rStrideD,1)/5), 1);
lStrideDExtracted = zeros(floor(size(lStrideD,1)/5), 1);

% re-initialize the iterators to 1
rj0 = 1;
lj0 = 1;

% find where the steps start in rStrideD
for ri = 1:size(smoothAccelR,1)
    if rV(ri,1) == 0 && rV(ri,2) == 0 && rV(ri,3) == 0
        rj0 = ri; % rj will be index of last set of zeros before step
    else
        break;
    end
end

% find where the steps start in lStrideD
for li = 1:size(smoothAccelL,1)
    if lV(li,1) == 0 && lV(li,2) == 0 && lV(li,3) == 0
        lj0 = li; % lj will be index of last set of zeros before step
    else
        break;
    end
end

% removing duplicates
rj = 2;
lj = 2;
rSize = 0;
lSize = 0;

for ri = rj0:size(rStrideD)-5
    if (rStrideD(ri,1)~=rStrideD(ri+1,1) && rStrideD(ri,1)~=rStrideD(ri+2,1) && rStrideD(ri,1)~=rStrideD(ri+3,1) && rStrideD(ri,1)~=rStrideD(ri+4,1) && rStrideD(ri,1)~=rStrideDExtracted(rj-1))
        rStrideDExtracted(rj) = rStrideD(ri,1);
        rj = rj + 1;
    end
    if rStrideD(ri,1) ~= 0
        rSize = rj;
    end
end

for li = lj0:size(lStrideD)-5
    if (lStrideD(li,1)~=lStrideD(li+1,1) && lStrideD(li,1)~=lStrideD(li+2,1) && lStrideD(li,1)~=lStrideD(li+3,1) && lStrideD(li,1)~=lStrideD(li+4,1) && lStrideD(li,1)~=lStrideDExtracted(lj-1))
        lStrideDExtracted(lj) = lStrideD(li,1);
        lj = lj + 1;
    end
    if lStrideD(li,1) ~= 0
        lSize = lj;
    end
end



% removing close duplicates
rStrideUndup = zeros(rSize,1);
lStrideUndup = zeros(lSize,1);
rj = 2;
lj = 2;

%Get individual right foot strides
for ri = 1:rSize
    if ~(abs(rStrideDExtracted(ri)) < abs(rStrideUndup(rj-1))+0.04 && abs(rStrideDExtracted(ri)) > abs(rStrideUndup(rj-1))-0.04)
        rStrideUndup(rj) = rStrideDExtracted(ri);
        rj = rj + 1;
    end
end

%Get individual left foot strides
for li = 1:lSize
    if ~(abs(lStrideDExtracted(li)) < abs(lStrideUndup(lj-1))+0.02 && abs(lStrideDExtracted(li)) > abs(lStrideUndup(lj-1))-0.02)
        lStrideUndup(lj) = lStrideDExtracted(li);
        lj = lj + 1;
    end
end

% remove 0s at end of array
for x = size(rStrideUndup):-1:2
    if rStrideUndup(x)==0
        rStrideUndup=rStrideUndup(1:end-1);
    end
end

% remove 0s at end of array
for x = size(lStrideUndup):-1:2
    if lStrideUndup(x)==0
        lStrideUndup=lStrideUndup(1:end-1);
    end
end

%use overall test length to increase accuracy of step measurements
measuredTotalR = rD(end,1);
measuredTotalL = abs(lD(end,1));

disp('Measured Total Right (m)');
disp(measuredTotalR);
disp('Measured Total Left (m)');
disp(measuredTotalL);

% get individual stride lengths from overall displacement
strideR = zeros(size(rStrideUndup));
for ri = 2:size(rStrideUndup)
    strideR(ri) = rStrideUndup(ri)-rStrideUndup(ri-1);
end

strideL = zeros(size(lStrideUndup));
for li = 2:size(lStrideUndup)
    strideL(li) = abs(lStrideUndup(li))-abs(lStrideUndup(li-1));
end

% use overall test length to increase accuracy of step measurements
correctedStrideR = (strideR/measuredTotalR)*len;
correctedStrideL = (strideL/measuredTotalL)*len;

% display stride data corrected using overall walk length
disp('---------------------');
disp('Corrected Right Foot Strides (m):');
disp(correctedStrideR);

disp('---------------------');
disp('Corrected Left Foot Strides (m):');
disp(correctedStrideL);
disp('---------------------');

%Get cumulative distance for each stride
correctedStrideRCum = zeros(size(correctedStrideR));
correctedStrideLCum = zeros(size(correctedStrideL));

for i = 2:size(correctedStrideR)
    correctedStrideRCum(i) = correctedStrideRCum(i-1)+correctedStrideR(i);
end

for i = 2:size(correctedStrideL)
    correctedStrideLCum(i) = correctedStrideLCum(i-1)+correctedStrideL(i);
end

% determining which foot stepped first
rightFirst = -1;
if correctedStrideLCum(2) < correctedStrideRCum(2) % left foot stepped first
    rightFirst = 0;
else % right foot stepped first
    rightFirst = 1;
end

% determining the largest size of array that we can use
if size(correctedStrideRCum,1) > size(correctedStrideLCum,1)
    max = size(correctedStrideLCum,1);
    stepLengthR = zeros(size(correctedStrideLCum));
    stepLengthL = zeros(size(correctedStrideLCum));
else
    max = size(correctedStrideRCum,1);
    stepLengthR = zeros(size(correctedStrideRCum));
    stepLengthL = zeros(size(correctedStrideRCum));
end

% getting the step lengths
if rightFirst == 1 % if stepped with right foot first
    for i = 2:max
        stepLengthR(i) = correctedStrideRCum(i) - correctedStrideLCum(i-1);
        stepLengthL(i) = correctedStrideLCum(i) - correctedStrideRCum(i);
    end
else % left foot stepped first
    for i = 2:max
        stepLengthL(i) = correctedStrideLCum(i) - correctedStrideRCum(i-1);
        stepLengthR(i) = correctedStrideRCum(i) - correctedStrideLCum(i);
    end
end

disp('Corrected Right Step Lengths (m):');
disp(stepLengthR);

disp('---------------------');
disp('Corrected Left Step Lengths (m):');
disp(stepLengthL);

%Plot X vs Y displacements
figure(8);
plot(rD(:,1),rD(:,2),'r',lD(:,1),lD(:,2),'b');
xlabel('Horizontal Displacement (m)');
ylabel('Vertical Displacement (m)');
title('Left And Right Foot Displacement');
