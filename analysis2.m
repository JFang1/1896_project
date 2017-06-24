%fileRight = fopen('right.txt','r');
%fileLeft = fopen('left.txt','r');

filename = 'Right.txt';
[accelR,delimeterOut]=importdata('Right.txt');

T = .023;
V = zeros(size(accelR,1),size(accelR,2));
D = zeros(size(accelR,1),size(accelR,2));


mag_accelR=abs(accelR);

Heel_Strike= mag_accelR< .1;


%Solve for velocity
% for y = 2:size(accelR,1)
%     for w = 1:size(accelR,2)
%         V(y,w) = V(y-1,w) + accelR(y,w)*T;     
%     end
% end

vel = zeros(size(accelR));
for t = 2:length(V)
    V(t,:) = V(t-1,:) + accelR(t,:) * T;
     if(Heel_Strike(t) == 1)
         V(t,:) = [0 0 0];     % force zero velocity when foot stationary
     end
end


% disp(V);

%Solve for displacement
for z = 2:size(accelR,1)
    for y = 1:size(accelR,2)
        D(z,y) = D(z-1,y) + V(z,y)*T;
    end
end

 disp(D);
%disp(mag_accelR);

plot3(D(:,1),D(:,2),D(:,3),'r')

% disp(sizeA);
% %Get the magnitude of the acceleration vector
% function f1 = mag(x,y,z)
%     f1 = sqrt(x^2 + y^2 + z^2);
% end

