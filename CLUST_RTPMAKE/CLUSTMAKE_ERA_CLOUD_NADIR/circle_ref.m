clear all
R = 10;

iMotion = input(' inertial frame (1) straight line or (2) SHM thru origin or (3) SHM in straight line away from origin : ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iMotion == 1
  t = 0 : 0.25 : 160;

  %% straight line motion, x constantly increasing, y = fixed
  x0 = -R;
  y0 = 0;

  x0 = +R;
  y0 = +R;

  v = -2/4;

  x_inertial = x0 + v*t;
  y_inertial = y0 * ones(size(t));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif iMotion == 2
  %% straight line motion through origin, x SHM, half period

  t = 0 : 0.1 : 20;
  x0 = +R;
  y0 = 0;
  x_inertial = x0 * cos(2*pi/40*t);
  y_inertial = y0 * ones(size(t));

elseif iMotion == 3
  %% straight line motion away from (0,0), x SHM, half period

  t = 0 : 0.1 : 20;
  x0 = +R;
  y0 = -R;
  x_inertial = x0 * cos(2*pi/40*t);
  y_inertial = y0 * ones(size(t));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf
scatter(x_inertial,y_inertial,30,t,'filled'); colormap default; hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = 40;
omega = 2*pi/T; 
for ii = 1 : length(t)
  matr = [cos(omega*t(ii)) +sin(omega*t(ii)); -sin(omega*t(ii)) cos(omega*t(ii))];
  inertial(1) = x_inertial(ii);
  inertial(2) = y_inertial(ii);  
  
  X = matr*inertial';
  x_rot(ii) = X(1);
  y_rot(ii) = X(2);
end

scatter(x_rot,y_rot,30,t,'filled'); hold off

colorbar; grid; axis equal; colormap jet