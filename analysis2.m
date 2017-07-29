clc
% NOTES:
%  - We MUST keep track of which file belongs to which foot
%  - We are assuming x and z are the dimensions of the floor/
%    horizontal plane (and y is the vertical dimension)
%  - x will be the direction of walking
%  - The patient MUST start with both feet at the same displacement

% request overall test length from user
len = input('Enter length of test (in meters): ');

% importing data for right foot
rightFile = 'RIGHT_POSVM3.txt';
[rAccel,rDelimeterOut] = importdata(rightFile);

% importing data for left foot
leftFile = 'LEFT_POSVM3.txt';
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
figure(1);
plot(lAccel(:,4),lAccel(:,1));
title('original left accel');

% define new array of X Y Z accel
ravg_x = rAccel(:,1);
ravg_y = rAccel(:,2);
ravg_z = rAccel(:,3);

lavg_x = lAccel(:,1);
lavg_y = lAccel(:,2);
lavg_z = lAccel(:,3);


stopsR = zeros(100,1);
stopsRCount=1;
stopsL = zeros(100,1);
stopsLCount=1;
startsR = zeros(100,1);
startsRCount=1;
startsL = zeros(100,1);
startsLCount=1;


% % filter X-axis acceleration vectors
% load('Lpass_Acc_X.mat');
% ravg_x = filter(NumX,1,rAccel(:,1));
% lavg_x = filter(NumX,1,lAccel(:,1));
%
% % filter Y-axis acceleration vectors
% load('Lpass_Acc_Y.mat');
% ravg_y = filter(NumAy,1,rAccel(:,2));
% lavg_y = filter(NumAy,1,lAccel(:,2));

% % plot filtered accel
% figure(2);
% plot(lAccel(:,4),lavg_x);
% title('filtered left accel')

% creating a smoothed acceleration curve using rolling average
movAvg = 14;
coeff = ones(1,movAvg)/movAvg;

% % create frequency spectrum of displacement
% Fs = 90;            % sampling frequency
% T = 1/Fs;             % sampling period
% L = length(lAccel(:,2));             % length of signal
% t = (0:L-1)*T;        % time vector
% lyf = fft(lAccel(:,2));
% P2 = abs(lyf/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = Fs*(0:(L/2))/L;
% figure(11);
% plot(f,P1);
% title('spectrum of displacement');

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

% plot filtered and smoothed curve
figure(3);
plot(lAccel(:,4),avg_x);
title('smoothed left signal');

% % find the heelstrikes
% smoothAccelR2 = smoothAccelR(:,1);
% smoothAccelL2 = smoothAccelL(:,1);
% [rpks, rlocs] = findpeaks(-smoothAccelR2);
% [lpks, llocs] = findpeaks(-smoothAccelL2);
% rHeelStrikes = zeros(size(smoothAccelR,1));
% lHeelStrikes = zeros(size(smoothAccelL,1));
% for i = 1:size(rlocs)
%     if (smoothAccelR2(i) < -0.4)
%         rHeelStrikes(rlocs(i)) = 1;
%     end
% end
% for i = 1:size(llocs)
%     if (smoothAccelL2(i) < -0.4)
%         lHeelStrikes(llocs(i)) = 1;
%     end
% end

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
        %detect if beginning of stride
        if(prevState==0)
            startsR(startsRCount)=w;
            startsRCount=startsRCount+1;
        end
        rV(w,:) = rV(w-1,:) + smoothAccelR(w,1:3) * rT(w)*(.001);
        % zero Out any backwards velocities
        if( rV(w,1) < 0)
            rV(w,:)=[0,0,0];
        end
        prevState = 1;
    end
    if(rHeelStrikes(w-2) == 1 && rHeelStrikes(w+2) == 1 && rHeelStrikes(w-1) == 1 && rHeelStrikes(w) == 1 && rHeelStrikes(w+1) == 1 || prevState==0)
        %Detect end of stride
        if(prevState==1)
            stopsR(stopsRCount)=w;
            stopsRCount=stopsRCount+1;
        end
        rV(w,:) = [0 0 0]; % force zero velocity when foot stationary
        prevState = 0;
    end
end

AVG_rV = mean(rV);

% calculating velocity data for the left foot
lAccelMag = abs(smoothAccelL);
lHeelStrikes = lAccelMag < .8;
lHeelLift = lAccelMag < .6;
prevState=0;
for w = 3:length(lV)-2
    % get time displacement for right foot
    lT(w)= smoothAccelL(w,4) - smoothAccelL(w-1,4);
    if(lHeelLift(w-2) == 0 && lHeelLift(w+2) == 0 && lHeelLift(w-1) == 0 && lHeelLift(w) == 0 && lHeelLift(w+1) == 0 || prevState==1)
        %detect if beginning of stride
        if(prevState==0)
            startsL(startsLCount)=w;
            startsLCount=startsLCount+1;
        end
        lV(w,:) = lV(w-1,:) + smoothAccelL(w,1:3) * lT(w)*(.001);
        % zero Out any backwards velocities
         if( lV(w,1) < 0)
             lV(w,:)=[0,0,0];
         end
        prevState=1;
    end
    if(lHeelStrikes(w-2) == 1 && lHeelStrikes(w+2) == 1 && lHeelStrikes(w-1) == 1 && lHeelStrikes(w) == 1 && lHeelStrikes(w+1) == 1 || prevState==0)
        %detect if end of stride
        if(prevState==1)
            stopsL(stopsLCount)=w;
            stopsLCount=stopsLCount+1;
        end
        lV(w,:) = [0 0 0]; % force zero velocity when foot stationary
        prevState = 0;
    end
end

AVG_lV = mean(lV);

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
%  2nd dimension (for z) is not needed if using only x for displacement
%  --keeping the 2nd dimension for now just in case we want to examine
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
        % rStrideD(rj,2) = rD(ri,3); % z dimension of displacement
        rj = rj + 1; % increment row in rStrideD
    end
end
% finding the displacements at heel strike - left foot
for li = 1:size(lStrideD,1) % for every sample
    if lV(li,1) == 0 && lV(li,2) == 0 && lV(li,3) == 0 % if heelstrike
        lStrideD(lj,1) = lD(li,1); % x dimension of displacement
        % lStrideD(lj,2) = lD(li,3); % z dimension of displacement
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

for ri = 1:rSize
    if ~(abs(rStrideDExtracted(ri)) < abs(rStrideUndup(rj-1))+0.04 && abs(rStrideDExtracted(ri)) > abs(rStrideUndup(rj-1))-0.04)
        rStrideUndup(rj) = rStrideDExtracted(ri);
        rj = rj + 1;
        % fprintf('O: %d %d\n', rStrideDExtracted(ri), rStrideUndup(rj-1));
    else
        % fprintf('X: %d %d\n', rStrideDExtracted(ri), rStrideUndup(rj-1));
    end
end

for li = 1:lSize
    if ~(abs(lStrideDExtracted(li)) < abs(lStrideUndup(lj-1))+0.02 && abs(lStrideDExtracted(li)) > abs(lStrideUndup(lj-1))-0.02)
        lStrideUndup(lj) = lStrideDExtracted(li);
        lj = lj + 1;
        % fprintf('O: %d %d\n', lStrideDExtracted(li), lStrideUndup(lj-1));
    else
        % fprintf('X: %d %d\n', lStrideDExtracted(li), lStrideUndup(lj-1));
    end
end

% remove 0s at end of array
for x = size(rStrideUndup):-1:2
    if rStrideUndup(x)==0
        rStrideUndup=rStrideUndup(1:end-1);
    end
end

for x = size(lStrideUndup):-1:2
    if lStrideUndup(x)==0
        lStrideUndup=lStrideUndup(1:end-1);
    end
end

% use overall test length to increase accuracy of step measurements
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

%Get number of strides for each foot
num_strR = find(~startsR);
num_strL = find(~startsL);

num_strR=num_strR(1);
num_strL=num_strL(1);

timeR=zeros(num_strR,1);
timeL=zeros(num_strL,1);

% display stride data corrected using overall walk length
disp('---------------------');
disp('Corrected Right Foot Strides:');
disp(sprintf('1: % 0.4fm',correctedStrideR(1)));
for str= 1:num_strR-1
   timeR(str) = smoothAccelR(stopsR(str),4)-smoothAccelR(startsR(str),4);
   if str<9
        newStr=sprintf('%d: % 0.4fm   %dms',str+1,correctedStrideR(str+1),timeR(str)); 
   else
        newStr=sprintf('%d: %0.4fm   %dms',str+1,correctedStrideR(str+1),timeR(str)); 
   end
   disp(newStr);
end

disp('---------------------');
disp('Corrected Left Foot Strides:');
disp(sprintf('1: % 0.4fm',correctedStrideL(1)));
for str= 1:num_strL-1
   timeL(str) = smoothAccelL(stopsL(str),4)-smoothAccelL(startsL(str),4);
   if str<9
       newStr=sprintf('%d: % 0.4fm   %dms',str+1,correctedStrideL(str+1),timeL(str)); 
   else
       newStr=sprintf('%d: %0.4fm   %dms',str+1,correctedStrideL(str+1),timeL(str)); 
   end
   disp(newStr);
end
disp('---------------------');

correctedStrideRCum = zeros(size(correctedStrideR));
correctedStrideLCum = zeros(size(correctedStrideL));

for i = 2:size(correctedStrideR)
    correctedStrideRCum(i) = correctedStrideRCum(i-1)+correctedStrideR(i);
end

for i = 2:size(correctedStrideL)
    correctedStrideLCum(i) = correctedStrideLCum(i-1)+correctedStrideL(i);
end

% disp(correctedStrideRCum);
% disp(correctedStrideLCum);

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
% 
% disp('max');
% disp(max);

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

disp('Right Step Length:');
for x = 1:max
    if x<10
        newStr=sprintf('%d: % 0.4fm',x,stepLengthR(x,1));
    else
        newStr=sprintf('%d: %0.4fm',x,stepLengthR(x,1));
    end
   disp(newStr);
end

disp('---------------------');
disp('Left Step Length:');

for x = 2:max
    if x<10
        newStr=sprintf('%d: % 0.4fm',x,stepLengthL(x,1));
    else
        newStr=sprintf('%d: %0.4fm',x,stepLengthL(x,1));
    end
   disp(newStr);
end


% plotting the step traces for both feet on the same graph
% figure(1);
% view(3);
% plot3(rD(:,1),rD(:,2),rD(:,3),'r');
% hold on;
% plot3(-lD(:,1),lD(:,2),lD(:,3),'b'); % if the data is the same, only the latter curve will appear
%
figure(8);
plot(rD(:,1),rD(:,2),'r',lD(:,1),lD(:,2),'b');
