% importing data for right foot
rightFile = 'Step1.TXT';
[rAccel,rDelimeterOut] = importdata(rightFile);

% importing data for left foot
leftFile = 'Step1.txt';
[lAccel,lDelimeterOut] = importdata(leftFile);

% constant
T = .023;

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
         rV(rt,:) = [0 0 0];     % force zero velocity when foot stationary
     end
end

% calculating velocity data for the left foot
lAccelMag = abs(lAccel);
lHeelStrikes = lAccelMag < .1;
lVel = zeros(size(lAccel));
for lt = 2:length(lV)
    lV(lt,:) = lV(lt-1,:) + lAccel(lt,:) * T;
     if(lHeelStrikes(lt) == 1)
         lV(lt,:) = [0 0 0];     % force zero velocity when foot stationary
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

% displaying displacement data
disp('------------------');
disp('------------------');
disp('Right Foot Displacement');
disp(rD);
disp('------------------');
disp('Left Foot Displacement');
disp(lD);

% displaying just the heel strikes
%TODO: find out why some of the heel strikes are 1
disp('------------------');
disp('Right Heel Strikes');
disp(rHeelStrikes);
disp('Left Heel Strikes');
disp(lHeelStrikes);

% plotting the steps for both feet on the same graph
%TODO: fix these so that they are on the same graph
plot3(rD(:,1),rD(:,2),rD(:,3),'r');
%plot3(lD(:,1),lD(:,2),lD(:,3),'b');




% 2nd method of finding where the heel strikes the ground (local minima)

% large number for initializing the size of the minima arrays
%maxSize = 2000;

% local minima of the Z-axis/vertical dimension, indicating heel strikes
%left_heel_strikes = zeros(max_size,3);
%right_heel_strikes = zeros(max_size,3);

% actually finding the minima
%left_heel_strikes = findpeaks(-leftData(:,3));
%right_heel_strikes = findpeaks(-rightData(:,3));