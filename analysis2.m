% NOTES: 
%  - We MUST keep track of which file belongs to which foot
%  - We are assuming x and y are the dimensions of the floor/
%    horizontal plane (and z is the vertical dimension)

% importing data for right foot
rightFile = 'RT_FOOT_21FT.TXT';
[rAccel,rDelimeterOut] = importdata(rightFile);

% importing data for left foot
leftFile = 'LT_FOOT_21FT.txt';
[lAccel,lDelimeterOut] = importdata(leftFile);

% constant
T = .012;

% initializing arrays that will hold velocity and displacement data
%  for the right foot
rV = zeros(size(rAccel,1),size(rAccel,2));
rD = zeros(size(rAccel,1),size(rAccel,2));
%  and for the left foot
lV = zeros(size(lAccel,1),size(lAccel,2));
lD = zeros(size(lAccel,1),size(lAccel,2));

% calculating velocity data for the right foot
rAccelMag = abs(rAccel);
rHeelStrikes = rAccelMag < .1;
rVel = zeros(size(rAccel));
for rt = 2:length(rV)
    rV(rt,:) = rV(rt-1,:) + rAccel(rt,:) * T;
     if(rHeelStrikes(rt) == 1)
         rV(rt,:) = [0 0 0 0 0];     % force zero velocity when foot stationary
     end
end

% calculating velocity data for the left foot
lAccelMag = abs(lAccel);
lHeelStrikes = lAccelMag < .1;
lVel = zeros(size(lAccel));
for lt = 2:length(lV)
    lV(lt,:) = lV(lt-1,:) + lAccel(lt,:) * T;
     if(lHeelStrikes(lt) == 1)
         lV(lt,:) = [0 0 0 0 0];     % force zero velocity when foot stationary
     end
end

% calculating displacement for the right foot
for rz = 2:size(rAccel,1)
    for ry = 1:size(rAccel,2)
        rD(rz,ry) = rD(rz-1,ry) + rV(rz,ry) * T;
    end
end

% calculating displacement for the left foot
for lz = 2:size(lAccel,1)
    for ly = 1:size(lAccel,2)
        lD(lz,ly) = lD(lz-1,ly) + lV(lz,ly) * T;
    end
end

% determining individual stride lengths . . . 

% array marker for right and left heel strike data in x+y dimensions
rStrideD = zeros(size(rAccel,1), 2);
lStrideD = zeros(size(lAccel,1), 2);
% separate iterators
rj = 1;
lj = 1;

% finding the displacements at heel strike - right foot
for ri = 2:size(rAccel,1) % for every sample
    if rV(ri,1) == 0 & rV(ri,2) == 0 & rV(ri,3) == 0 % if heelstrike
        rStrideD(rj,1) = rD(ri,1); % x dimension of displacement
        rStrideD(rj,2) = rD(ri,2); % y dimension of displacement
        rj = rj + 1; % increment row in rStrideD
    end
end
% finding the displacements at heel strike - left foot
for li = 2:size(lAccel,1) % for every sample
    if lV(li,1) == 0 & lV(li,2) == 0 & lV(li,3) == 0 % if heelstrike
        lStrideD(lj,1) = lD(li,1); % x dimension of displacement
        lStrideD(lj,2) = lD(li,2); % y dimension of displacement
        lj = lj + 1; % increment row in lStrideD
    end
end

% determining individual step lengths . . .

% re-initialize the iterators to 1
rj = 1;
lj = 1;

% create vectors for step distance

% displaying displacement data
disp('------------------');
disp('------------------');
%disp('Right Foot Displacement');
disp(rStrideD);
disp('------------------');
disp('Left Foot Displacement');
disp(lStrideD);

% plotting the step traces for both feet on the same graph
figure(1);
view(3);
plot3(rD(:,1),rD(:,2),rD(:,3),'r');
hold on;
plot3(-lD(:,1),lD(:,2),lD(:,3),'b'); % if the data is the same, only the latter curve will appear




% 2nd method of finding where the heel strikes the ground (local minima)

% large number for initializing the size of the minima arrays
%maxSize = 2000;

% local minima of the Z-axis/vertical dimension, indicating heel strikes
%left_heel_strikes = zeros(max_size,3);
%right_heel_strikes = zeros(max_size,3);

% actually finding the minima
%left_heel_strikes = findpeaks(-leftData(:,3));
%right_heel_strikes = findpeaks(-rightData(:,3));
figure(2);
plot(rD(:,1),rD(:,2),'r',-lD(:,1),lD(:,2),'b');





