%fileRight = fopen('right.txt','r');
%fileLeft = fopen('left.txt','r');

filename = 'Right.txt';
[accelR,delimeterOut]=importdata('Right.txt');

T = .5;
V = zeros(size(accelR,1),size(accelR,2));
D = zeros(size(accelR,1),size(accelR,2));

%Solve for velocity
for y = 1:size(accelR,1)
    for w = 1:size(accelR,2)
        V(y,w) = 0 + accelR(y,w)*T;
    end
end

disp(V);

%Solve for displacement
for z = 1:size(accelR,1)
    for y = 1:size(accelR,2)
        D(z,y) = 0 + V(z,y)*T;
    end
end

disp(D);
disp(sizeA);
%Get the magnitude of the acceleration vector
function f1 = mag(x,y,z)
    f1 = sqrt(x^2 + y^2 + z^2);
end

