% NOTES: 
%  - We MUST keep track of which file belongs to which foot
%  - We are assuming x and z are the dimensions of the floor/
%    horizontal plane (and y is the vertical dimension)
%  - x will be the direction of walking
%  - The patient MUST start with both feet at the same displacement

%Request overall test length from user
len = input('Enter Length of Test (in meters): ');

% importing data for right foot
rightFile = 'RT_FOOT_21FT.txt';
[rAccel,rDelimeterOut] = importdata(rightFile);

% importing data for left foot
leftFile = 'LT_FOOT_21FT.txt';
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

%Creating a smoothed Acceleration curve using rolling average

movAvg = 6;
coeff = ones(1,movAvg)/movAvg;

%Smooth Right foot curve
avg_x = filter(coeff,1,rAccel(:,1));
avg_y = filter(coeff,1,rAccel(:,2));
avg_z = filter(coeff,1,rAccel(:,3));
smoothAccelR = [avg_x,avg_y,avg_z, rAccel(:,4)];

%Smooth Left foot curve
avg_x = filter(coeff,1,lAccel(:,1));
avg_y = filter(coeff,1,lAccel(:,2));
avg_z = filter(coeff,1,lAccel(:,3));
smoothAccelL = [avg_x,avg_y,avg_z, lAccel(:,4)];


% calculating velocity data for the right foot
rAccelMag = abs(smoothAccelR);
rHeelStrikes = rAccelMag(:,1) < .4;
for w = 2:length(rV)-1
    %get time displacement for right foot
    rT(w)= smoothAccelR(w,4) - smoothAccelR(w-1,4);
    rV(w,:) = rV(w-1,:) + smoothAccelR(w,1:3) * rT(w)*(.001);
     if(rHeelStrikes(w-1) == 1 && rHeelStrikes(w) == 1 && rHeelStrikes(w+1) == 1)
         rV(w,:) = [0 0 0];     % force zero velocity when foot stationary
     end
end

% calculating velocity data for the left foot
lAccelMag = abs(smoothAccelL);
lHeelStrikes = lAccelMag(:,1) < .4;
for w = 2:length(lV)-1
    %get time displacement for right foot
    lT(w)= smoothAccelL(w,4) - smoothAccelL(w-1,4);
    lV(w,:) = lV(w-1,:) + smoothAccelL(w,1:3) * lT(w)*(.001);
     if(lHeelStrikes(w-1) == 1 && lHeelStrikes(w) == 1 && lHeelStrikes(w+1) == 1)
         lV(w,:) = [0 0 0]; % force zero velocity when foot stationary
     end
end

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


% determining individual stride lengths . . . 

% array marker for right and left heel strike data in x+z dimensions
%  2nd dimension (for z) is not needed if using only x for displacement
%  I'm keeping the 2nd dimension for now just in case we want to examine
%  other calculation results.
rStrideD = zeros(size(rAccel,1), 2);
lStrideD = zeros(size(lAccel,1), 2);
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
for li = 2:size(lStrideD,1) % for every sample
    if lV(li,1) == 0 && lV(li,2) == 0 && lV(li,3) == 0 % if heelstrike
        lStrideD(lj,1) = lD(li,1); % x dimension of displacement
        % lStrideD(lj,2) = lD(li,3); % z dimension of displacement
        lj = lj + 1; % increment row in lStrideD
    end
end

% Getting rid of repetition in the stride data

rStrideDExtracted = zeros(floor(size(rStrideD,1)/5), 1);
lStrideDExtracted = zeros(floor(size(lStrideD,1)/5), 1);

% re-initialize the iterators to 1
rj0 = 1;
lj0 = 1;

% find where the steps start in rStrideD
for ri = 1:size(rAccel,1)
<<<<<<< HEAD
    if rV(ri,1) == 0 && rV(ri,2) == 0 && rV(ri,3) == 0
=======
    if rV(ri,1) == 0 & rV(ri,2) == 0 & rV(ri,3) == 0
>>>>>>> 33eba1c4130bb2a831596314d980f96a73811ead
        rj0 = ri; % rj will be index of last set of zeros before step
    else
        break;
    end
end

% find where the steps start in lStrideD
for li = 1:size(lAccel,1)
<<<<<<< HEAD
    if lV(li,1) == 0 && lV(li,2) == 0 && lV(li,3) == 0
=======
    if lV(li,1) == 0 & lV(li,2) == 0 & lV(li,3) == 0
>>>>>>> 33eba1c4130bb2a831596314d980f96a73811ead
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

for ri = rj0+5:size(rStrideD)
<<<<<<< HEAD
    if (rStrideD(ri,1)~=rStrideD(ri-1,1) && rStrideD(ri,1)~=rStrideD(ri-2,1) && rStrideD(ri,1)~=rStrideD(ri-3,1) && rStrideD(ri,1)~=rStrideD(ri-4,1) && rStrideD(ri,1)~=rStrideDExtracted(rj-1))
=======
    if (rStrideD(ri,1)~=rStrideD(ri-1,1) & rStrideD(ri,1)~=rStrideD(ri-2,1) & rStrideD(ri,1)~=rStrideD(ri-3,1) & rStrideD(ri,1)~=rStrideD(ri-4,1) & rStrideD(ri,1)~=rStrideDExtracted(rj-1))
>>>>>>> 33eba1c4130bb2a831596314d980f96a73811ead
        rStrideDExtracted(rj) = rStrideD(ri,1);
        rj = rj + 1;
    end
    if rStrideD(ri,1) ~= 0
        rSize = rj;
    end
end
    
for li = lj0+5:size(lStrideD)
<<<<<<< HEAD
    if (lStrideD(li,1)~=lStrideD(li-1,1) && lStrideD(li,1)~=lStrideD(li-2,1) && lStrideD(li,1)~=lStrideD(li-3,1) && lStrideD(li,1)~=lStrideD(li-4,1) && lStrideD(li,1)~=lStrideDExtracted(lj-1))
=======
    if (lStrideD(li,1)~=lStrideD(li-1,1) & lStrideD(li,1)~=lStrideD(li-2,1) & lStrideD(li,1)~=lStrideD(li-3,1) & lStrideD(li,1)~=lStrideD(li-4,1) & lStrideD(li,1)~=lStrideDExtracted(lj-1))
>>>>>>> 33eba1c4130bb2a831596314d980f96a73811ead
        lStrideDExtracted(lj) = lStrideD(li,1);
        lj = lj + 1;
    end
    if lStrideD(li,1) ~= 0
        lSize = lj;
    end
end

% removing close duplicates
rStrideUndup = zeros(floor(rSize/5),1);
lStrideUndup = zeros(floor(lSize/5),1);
rj = 2;
lj = 2;

for ri = 2:rSize
<<<<<<< HEAD
    if ~(abs(rStrideDExtracted(ri)) < abs(rStrideUndup(rj-1))+0.04 && abs(rStrideDExtracted(ri)) > abs(rStrideUndup(rj-1))-0.04)
=======
    if ~(abs(rStrideDExtracted(ri)) < abs(rStrideUndup(rj-1))+0.04 & abs(rStrideDExtracted(ri)) > abs(rStrideUndup(rj-1))-0.04)
>>>>>>> 33eba1c4130bb2a831596314d980f96a73811ead
        rStrideUndup(rj) = rStrideDExtracted(ri);
        rj = rj + 1;
        fprintf('O: %d %d\n', rStrideDExtracted(ri), rStrideUndup(rj-1));
    else
        fprintf('X: %d %d\n', rStrideDExtracted(ri), rStrideUndup(rj-1));
    end
end
    
for li = 2:lSize
<<<<<<< HEAD
    if ~(abs(lStrideDExtracted(li)) < abs(lStrideUndup(lj-1))+0.02 && abs(lStrideDExtracted(li)) > abs(lStrideUndup(lj-1))-0.02)
=======
    if ~(abs(lStrideDExtracted(li)) < abs(lStrideUndup(lj-1))+0.02 & abs(lStrideDExtracted(li)) > abs(lStrideUndup(lj-1))-0.02)
>>>>>>> 33eba1c4130bb2a831596314d980f96a73811ead
        lStrideUndup(lj) = lStrideDExtracted(li);
        lj = lj + 1;
        fprintf('O: %d %d\n', lStrideDExtracted(li), lStrideUndup(lj-1));
    else
        fprintf('X: %d %d\n', lStrideDExtracted(li), lStrideUndup(lj-1));
    end
end
<<<<<<< HEAD

%remove 0 at end of array
rStrideUndup = rStrideUndup(1,1:end-1);
lStrideUndup = abs(lStrideUndup(1,1:end-1));
=======
>>>>>>> 33eba1c4130bb2a831596314d980f96a73811ead

disp('---------------------');
disp('---------------------');
disp('Right foot heelstrike displacements without repetition:');
<<<<<<< HEAD
disp(rStrideUndup.');
disp('---------------------');
disp('Left foot heelstrike displacements without repetition:');
disp(lStrideUndup.');
disp('---------------------');

%Use overall test length to increase accuracy of step measurements
measuredTotalR = rD(end,1);
measuredTotalL = abs(lD(end,1));

%Get individual stride lengths from overall displacement
strideR = zeros(size(rStrideUndup.'));
for ri = 2:size(rStrideUndup.')
    strideR(ri) = rStrideUndup(ri)-rStrideUndup(ri-1);
end

strideL = zeros(size(lStrideUndup.'));
for li = 2:size(lStrideUndup.')
    strideL(li) = lStrideUndup(li)-lStrideUndup(li-1);
end

%Use overall test length to increase accuracy of step measurements
correctedStrideR = (strideR/measuredTotalR)*len;
correctedStrideL = (strideL/measuredTotalL)*len;

%Display stride data corrected using overall walk length
disp('---------------------');
disp('Corrected Right foot heelstrike displacements without repetition:');
disp(correctedStrideR);
disp('---------------------');
disp('Corrected Left foot heelstrike displacements without repetition:');
disp(correctedStrideL);
disp('---------------------');

=======
disp(rStrideUndup);
disp('---------------------');
disp('Left foot heelstrike displacements without repetition:');
disp(lStrideUndup);
disp('---------------------');

>>>>>>> 33eba1c4130bb2a831596314d980f96a73811ead
% determining individual step lengths . . .

% % converting cumulative distances to individual stride lengths along x,y
% for ri = size(rD,1)-rj:1
%     rD(ri,1) = rD(ri,1) - rD(ri-1,1);
%     rD(ri,2) = rD(ri,2) - rD(ri-1,2);
% end
% 
% for li = size(lD,1)-lj:1
%     lD(li,1) = lD(li,1) - lD(li-1,1);
%     lD(li,2) = lD(li,2) - lD(li-1,2);
% end
% 
% % create vectors for stride distance in 1 dimension
% % NOTE: This assumes walking in a straight line, no turning
% rStrideD = zeros(size(rAccel,1)-rj+1, 2); % mathematically correct to add 1
% lStrideD = zeros(size(lAccel,1)-lj+1);
% 
% for ri = 1:size(rStrideD)
%     rStrideD(ri) = sqrt(rD(rj,1)^2 + rD(rj,2)^2);
% end
% 
% for li = 1:size(lStrideD)
%     lStrideD(li) = sqrt(lD(lj,1)^2 + lD(lj,2)^2);
% end
% 
% % determininig which foot stepped first
% temp = 0;
% rightFirst = -1;
% if lStrideD(1) < rStrideD(1) % left foot stepped first
%     rightFirst = 0;
% else % right foot stepped first
%     rightFirst = 1;
% end
% 
% if rightFirst == 1 % if stepped with right foot first
%     for li = 2:2:size(lStrideD)
%         lStrideD(li) = lStrideD(li) - rStrideD(li);
%         rStrideD(li+1) = rStrideD(li+1) - lStrideD(li);
%     end
% else
%     
% end


% plotting the step traces for both feet on the same graph
% figure(1);
% view(3);
% plot3(rD(:,1),rD(:,2),rD(:,3),'r');
% hold on;
% plot3(-lD(:,1),lD(:,2),lD(:,3),'b'); % if the data is the same, only the latter curve will appear
% 
% figure(2);
<<<<<<< HEAD
plot(rD(:,1),rD(:,2),'r',-lD(:,1),lD(:,2),'b');
=======
% plot(rD(:,1),rD(:,2),'r',-lD(:,1),lD(:,2),'b');
>>>>>>> 33eba1c4130bb2a831596314d980f96a73811ead
