% Forced damping 1 DOF vibrating system
clear
clc
%% Given values
F0 = 50;
m = 5;
c = 10;
k = 250;
%% Initial condition
x0 = 0;
v0 = 0;
%% Calculation
z = c/(2*sqrt(m*k))
wn = sqrt(k/m);
r=sqrt(1-2*z^2) % Frequency ratios
time = linspace(0,4,200);

%% Simulation 1 : Displacement response
        wd = sqrt(1-z^2)*wn;
        M=1/sqrt((1-r^2)^2+(2*z*r)^2)
        X = F0/(k*M)
        w = r*wn;
        phi = atan((2*z*r)/(1-r^2))*180/pi
figure;
syms x(t);            % Define dependent variable
Dx = diff(x,t);                                   % Define first derivation
D2x = diff(Dx,t);                                  % Define second derivation                                    

x_res = dsolve( m*D2x + c*Dx + k*x - F0*sin(w*t) == 0  ,x(0)==x0 ,Dx(0) == v0,  't');  % Initial conditions units are [m] and [m/s] respectively
vx_res = diff(x_res);

x_fun = matlabFunction(x_res);                    % Create function handle for x i.e converting to matlab function
vx_fun = matlabFunction(vx_res);               % Create function handle for x_dot

x_t = x_fun(time);                                % Evaluate function at points "t"
vx_t = vx_fun(time);                              % Evaluate function at points "t"

%% Plotting
yyaxis left
plot(time, x_t,'LineWidth',1.5, 'MarkerSize', 10);
ylabel('Displacement (m)');
xlabel('Time(s)');
x_t_max_limit = (max(abs(x_t))+0.05*max(abs(x_t))); % Get Maximum of x_t and add 5 percent
ylim([-x_t_max_limit,x_t_max_limit]);
yyaxis right
plot(time,vx_t,'LineWidth',1.5, 'MarkerSize', 10);
ylabel('Velocity (m/s)');
vx_t_max_limit = (max(abs(vx_t))+0.05*max(abs(vx_t))); % Get Maximum of v_t and add 5 percent
ylim([-vx_t_max_limit,vx_t_max_limit]); 
title('Forced Vibration');
grid on


%% Simulation 2 : Magnification factor
figure;
r = linspace(0,2,100);
Z = linspace(0,1,5); 
for i=1:length(Z)
    M=1./sqrt((1-r.^2).^2+(2*Z(i)*r).^2);
    leg2{i}=sprintf('zeta = %0.2f',Z(i));
    subplot(2,2,1);
    plot(r,M,'LineWidth',1.5, 'MarkerSize', 10)   
    hold on 
end
legend(leg2);
title('Magnification plot (for various zeta)');
xlabel({'Frequency ratio = r'});
ylabel({'Magnification factor = M'});
ylim([0,4]);
grid on;
hold off;
subplot(2,2,2);
M0=1./sqrt((1-r.^2).^2+(2*z*r).^2);
plot(r,M0,'LineWidth',1.5, 'MarkerSize', 10);
hold on;
leg3=sprintf('zeta = %0.4f',z);
legend(leg3);
title('Magnification plot for zeta of problem');
xlabel({'Frequency ratio = r'});
ylabel({'Magnification factor = M'});
ylim([0,4]);
grid on;
hold off;


%% Simulation 3 : phi
for i=1:length(Z)
    phi_full = atan((2.*Z(i).*r)./(1-r.^2)).*180./pi;
    subplot(2,2,3);
    plot(r,phi_full,'LineWidth',1.5, 'MarkerSize', 10);
    hold on 
end
legend(leg2);
title('phi plot (for various zeta)');
xlabel({'Frequency ratio = r'});
ylabel({'phi','(deg)'});
ylim([0,180]);
grid on;
hold off;

theta = atan((2.*z.*r)./(1-r.^2)).*180./pi;
subplot(2,2,4);
plot(r,theta,'LineWidth',1.5, 'MarkerSize', 10);
hold on;
legend(leg3);
title('phi plot for zeta of problem');
xlabel({'Frequency ratio = r'});
ylabel({'phi','(deg)'});
ylim([0,180]);
xlim([0,2]);
grid on;
hold off;