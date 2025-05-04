"%% HomMat3D matrix copied from Halcon"
HomMat3D = [R];
H = reshape(HomMat3D', [], 1); % Change the matrix shape for convenience, no practical meaning.
%% Calculate the rotation axis vector and rotation angle
theta2u = H(10) - H(7);
theta2v = H(3) - H(9);
theta2w = H(5) - H(2)
u = 1 / sqrt(1 + (theta2v / theta2u)^2 + (theta2w / theta2u)^2);
v = u * theta2v / theta2u;
w = u * theta2w / theta2u;
sintheta = theta2u / (2 * u);
costheta = (H(11) - w^2) / (1 - w^2);
%% Compute the rotation axis matrix
A = [(1 - u^2) * (1 - costheta), w * sintheta - u * v * (1 - costheta), -v * sintheta - u * w * (1 - costheta);
     -w * sintheta - u * v * (1 - costheta), (1 - v^2) * (1 - costheta), u * sintheta - v * w * (1 - costheta);
     v * sintheta - v * w * (1 - costheta), -u * sintheta - v * w * (1 - costheta), (1 - w^2) * (1 - costheta)];
y = [H(4), H(8), H(12)];
As = A(:, 2:3);
%% Output result as [x, y, z]' = k * a' + b0' (k can take any value)
zz = [0; pinv(As) * y'];
ss = [u, v, w];
angle_rad = atan2(sintheta, costheta); % Compute angle in radians
angle_deg = rad2deg(angle_rad);       % Convert angle to degrees
angle_deg
%% Select two points on the rotation axis
x1 = 1;
y1 = (x1 - zz(1)) / ss(1) * ss(2) + zz(2);
z1 = (x1 - zz(1)) / ss(1) * ss(3) + zz(3);
x_1 = [x1, y1, z1];
x2 = 2;
y2 = (x2 - zz(1)) / ss(1) * ss(2) + zz(2);
z2 = (x2 - zz(1)) / ss(1) * ss(3) + zz(3);
x_2 = [x2, y2, z2];
%% Any point in space
x_3 = [R];
%% Calculate the foot of the perpendicular (N)
k = -((x_1(1) - x_3(1)) * (x_2(1) - x_1(1)) + ...
      (x_1(2) - x_3(2)) * (x_2(2) - x_1(2)) + ...
      (x_1(3) - x_3(3)) * (x_2(3) - x_1(3))) / (norm(x_1 - x_2)^2);
xn = k * (x_2(1) - x_1(1)) + x_1(1); xn = vpa(xn);
yn = k * (x_2(2) - x_1(2)) + x_1(2); yn = vpa(yn);
zn = k * (x_2(3) - x_1(3)) + x_1(3); zn = vpa(zn);
N = [xn, yn, zn];
fprintf("\nCoordinates of the foot of the perpendicular: %4.4f, %4.4f, %4.4f", xn, yn, zn);
% Distance formula
dis = sqrt((x_3(1) - N(1, 1))^2 + (x_3(2) - N(1, 2))^2 + (x_3(3) - N(1, 3))^2);
fprintf("\nDistance to the foot of the perpendicular: %4.4f", dis);
%% Compute the rotation direction vector
% Vector from N to X3
v_X3 = x_3 - N;
% Rotation direction vector (cross product)
v_rot = cross(ss, v_X3);
% Normalize the rotation direction vector
v_rot_norm = v_rot / norm(v_rot);
% Output the rotation direction vector
fprintf("\nRotation direction vector: %4.4f, %4.4f, %4.4f\n", v_rot_norm(1), v_ro