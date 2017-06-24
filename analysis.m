%fileRight = fopen('right.txt','r');
%fileLeft = fopen('left.txt','r');

filename = 'Right.txt';
[accelR,delimeterOut]=importdata('Right.txt');

T = .5;
V = zeros(size(accelR,1),size(accelR,2));
D = zeros(size(accelR,1),size(accelR,2));




%Solve for velocity
for y = 2:size(accelR,1)
    for w = 1:size(accelR,2)
        V(y,w) = V(y-1,w) + accelR(y,w)*T;     
    end
end

disp(V);

%Solve for displacement
for z = 2:size(accelR,1)
    for y = 1:size(accelR,2)
        D(z,y) = D(z-1,y) + V(z,y)*T;
    end
end

disp(D);

plot3(D(:,1),D(:,2),D(:,3),'r')

% disp(sizeA);
% %Get the magnitude of the acceleration vector
% function f1 = mag(x,y,z)
%     f1 = sqrt(x^2 + y^2 + z^2);
% end

