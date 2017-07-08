% NOTES: 
%  - We MUST keep track of which file belongs to which foot
%  - We are assuming x and z are the dimensions of the floor/
%    horizontal plane (and y is the vertical dimension)
%  - x will be the direction of walking
%  - The patient MUST start with both feet at the same displacement

% importing data for right foot
rightFile = 'RT_FOOT_21FT.TXT';
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

% calculating velocity data for the right foot
rAccelMag = abs(rAccel);
rHeelStrikes = rAccelMag(:,1) < .1;
for w = 2:length(rV)
    %get time displacement for right foot
    rT(w)= rAccel(w,4) - rAccel(w-1,4);
    rV(w,:) = rV(w-1,:) + rAccel(w,1:3) * rT(w)*(.001);
     if(rHeelStrikes(w) == 1)
         rV(w,:) = [0 0 0];     % force zero velocity when foot stationary
     end
end

% calculating velocity data for the left foot
lAccelMag = abs(lAccel);
lHeelStrikes = lAccelMag(:,1) < .1;
for w = 2:length(lV)
    %get time displacement for right foot
    lT(w)= rAccel(w,4) - rAccel(w-1,4);
    lV(w,:) = lV(w-1,:) + lAccel(w,1:3) * lT(w)*(.001);
     if(lHeelStrikes(w) == 1)
         lV(w,:) = [0 0 0]; % force zero velocity when foot stationary
     end
end

% calculating displacement for the right foot
for ri = 2:size(rAccel,1)
    for ry = 1:3
        rD(ri,ry) = rD(ri-1,ry) + rV(ri,ry) * rT(ri)*(.001);   
    end
end

% calculating displacement for the left foot
for li = 2:size(lAccel,1)
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
    if rV(ri,1) == 0 & rV(ri,2) == 0 & rV(ri,3) == 0 % if heelstrike
        rStrideD(rj,1) = rD(ri,1); % x dimension of displacement
        % rStrideD(rj,2) = rD(ri,3); % z dimension of displacement
        rj = rj + 1; % increment row in rStrideD
    end
end
% finding the displacements at heel strike - left foot
for li = 2:size(lStrideD,1) % for every sample
    if lV(li,1) == 0 & lV(li,2) == 0 & lV(li,3) == 0 % if heelstrike
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
    if rV(ri,1) == 0 & rV(ri,2) == 0 & rV(ri,3) == 0
        rj0 = ri; % rj will be index of last set of zeros before step
    else
        break;
    end
end

% find where the steps start in lStrideD
for li = 1:size(lAccel,1)
    if lV(li,1) == 0 & lV(li,2) == 0 & lV(li,3) == 0
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
    if (rStrideD(ri,1)~=rStrideD(ri-1,1) & rStrideD(ri,1)~=rStrideD(ri-2,1) & rStrideD(ri,1)~=rStrideD(ri-3,1) & rStrideD(ri,1)~=rStrideD(ri-4,1) & rStrideD(ri,1)~=rStrideDExtracted(rj-1))
        rStrideDExtracted(rj) = rStrideD(ri,1);
        rj = rj + 1;
    end
    if rStrideD(ri,1) ~= 0
        rSize = rj;
    end
end
    
for li = lj0+5:size(lStrideD)
    if (lStrideD(li,1)~=lStrideD(li-1,1) & lStrideD(li,1)~=lStrideD(li-2,1) & lStrideD(li,1)~=lStrideD(li-3,1) & lStrideD(li,1)~=lStrideD(li-4,1) & lStrideD(li,1)~=lStrideDExtracted(lj-1))
        lStrideDExtracted(lj) = lStrideD(li,1);
        lj = lj + 1;
    end
    if lStrideD(li,1) ~= 0
        lSize = lj;
    end
end

% removing close duplicates
rStrideUndup = zeros(floor(rSize/3),1);
lStrideUndup = zeros(floor(lSize/3),1);
rj = 2;
lj = 2;

for ri = 2:rSize
    if ~(abs(rStrideDExtracted(ri)) < abs(rStrideUndup(rj-1))+0.04 & abs(rStrideDExtracted(ri)) > abs(rStrideUndup(rj-1))-0.04)
        rStrideUndup(rj) = rStrideDExtracted(ri);
        rj = rj + 1;
        fprintf('O: %d %d\n', rStrideDExtracted(ri), rStrideUndup(rj-1));
    else
        fprintf('X: %d %d\n', rStrideDExtracted(ri), rStrideUndup(rj-1));
    end
end
    
for li = 2:lSize
    if ~(abs(rStrideDExtracted(li)) < abs(lStrideUndup(lj-1))+0.02 & abs(lStrideDExtracted(li)) > abs(lStrideUndup(lj-1))-0.02)
        lStrideUndup(lj) = lStrideDExtracted(li);
        lj = lj + 1;
        fprintf('O: %d %d\n', lStrideDExtracted(li), lStrideUndup(lj-1));
    else
        fprintf('X: %d %d\n', lStrideDExtracted(li), lStrideUndup(lj-1));
    end
end

disp('---------------------');
disp('---------------------');
disp('Right foot heelstrike displacements without repetition:');
disp(rStrideUndup);
disp('---------------------');
disp('Left foot heelstrike displacements without repetition:');
disp(lStrideUndup);
disp('---------------------');

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
figure(1);
view(3);
plot3(rD(:,1),rD(:,2),rD(:,3),'r');
hold on;
plot3(-lD(:,1),lD(:,2),lD(:,3),'b'); % if the data is the same, only the latter curve will appear

figure(2);
plot(rD(:,1),rD(:,2),'r',-lD(:,1),lD(:,2),'b');