% matlab 2022b
% author: zhaochichi github: https://github.com/chichizhao/clone_jzg22
% function: plot the 3D figure of the clone mutation pattern
%set x axis to 86
xlim([0 86])
%set y axis to 72
ylim([0 72])
t =[]
for i = 1:70 
    t(i,1) = -1400;
end
% add a 0 to x_points
x_points = [0;x_points];
% add a 0 to y_points
y_points = [0;y_points];
% add a 0 to z
z = [0;z];
line= plot3(y_points,x_points,t)
%set the linewidth to 2 
set(line,'linewidth',2)

x1_points = []
for i = 1:70 
    x1_points(i+i-1,1) = x_points(i);
    x1_points(i+i,1) = x_points(i);
end
y1_points = []
for i = 1:70 
    y1_points(i+i-1,1) = y_points(i);
    y1_points(i+i,1) = y_points(i);
end
z1_points = []
for i = 1:70 
    z1_points(i+i-1,1) = z(i);
    z1_points(i+i,1) = z(i);
end

for i = 1:length(x1_points)
    t(i,1) = 72;
end
line=plot3(y1_points+1,t,-z1_points+0.5)
%set the linewidth to 2
set(line,'linewidth',2)

for i = 1:length(y1_points)
    t(i,1) = 86;
end
line= plot3(t,x1_points+1,-z1_points+0.5)
%set the linewidth to 2
set(line,'linewidth',2)

line=plot3(y1_points+1,x1_points+1,-z1_points+3)
%set the linewidth to 2
set(line,'linewidth',2)

z_mesh2 = zeros(172, 172);
for i = 1:172
    for j = 1:172
        z_mesh2(i,j) = -1400;
    end
end
surf(x,y,z_mesh2,c);
%set the z axis to 2400
zlim([0 -2400])
%set the z axis to right
set(gca,'ZDir','reverse')
xlabel('mutation density')