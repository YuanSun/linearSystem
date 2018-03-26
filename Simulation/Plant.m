% plant.m
% This is the script Calculates inertias and builds the rigid body or
% flexible plant models
% Yuan Sun
% 2015-11-20

clc
clear

%% ****  Input Parameters ****

%%
% Select output location, Out:
%
%        % "1" = Encoder #1; 
%        % "2" = Encoder #2
Out = 2;    % choose Encoder #2

%%
% Plant configuration, cofig:
%
%        % "1" = Drive disk only, (use only out = 1)
%        % "2" = Drive & load disk -- rigid drive train
%        % "3" = Drive & load disk -- flexixble drive belt

Config = 2; % choose flexible drive belt

%%
% Selec mass parameters
mwd = 0.8; %(mass of brass weights on drive disk)
rwd = 0.05; %(radius from center of plate on drive disk)
mwl = 2.0; %(mass of brass weight on load disk)
rwl = 0.1; %(radius from center of plate on load disk)

%Lower (Drive) Pully: Select (decomment) one of the following as
%appropriate:
npd = 18; Jpd = 0.000003;  %(18 tooth)
%npd = 24, Jpd = 0.000008;   %(24 tooth)
%npd = 36, Jpd = 0.000039;  %(36 tooth)
%npd = 72, Jpd = 0.00055;   %(72 tooth)

%Upper (Load) Pulley: Select (decomment) one of the followint as
%appropriate:
%npl = 18, Jpl = 0.000003;  %(18 tooth)
%npl = 24, Jpl = 0.000008;  %(24 tooth)
%npl = 36, Jpl = 0.000039;  %(36 tooth)
npl = 72; Jpl = 0.00055;    %(72 tooth)

%Default frict coeff's -- change if setup-specific measurements available
c1 = 0.004;  % if belt from drive to SR ass'y attached
%c1 = 0.002 % if belt from drive to SR ass'y not attached
c2 = 0.05;   % if belt from drive to SR ass'y attached
c12 = 0.017; % viscous coupling between drive & load

k = 8.45; % Torsional spring constant

khw = 5.81;  %(Hardware gain, assume kg=1, it will be corrected below if out = 2)


%% ****  Transfer Function and State Space Model  ****

%%
% Gear ratios

gr = 6 * npd / npl; %24;  % This is the value we found    % 6 * npd / npl; --> this is suggested by instructor  % gear ratio
grprime = npd / 12; %6; % This is the value we found  % npd / 12; --> This is suggested % drive to SR pulley gear ratio

%%
% First calculate known inertias
Jdd = 0.0004;         %0.000326; % 0.000400 is given by instructor's suggestion
Jdl = 0.0065;
Jpbl = 0.000031;    % Backlash mechanism
Jp = Jpd+Jpl+Jpbl;

%%
% Calculate Drive inertia

% initializing
rwdo = 0; 
rwlo = 0;

%%
% Find which size weight used, can use only 0, 2, or 4 weights
if mwd < 0.81
    if mwd < 0.39
        rwdo = 0.016; % (smaller brass weight used)
    end
end

if mwd < 2.1
    if mwd > 0.9
        rwdo = 0.025;   % (larger brass weight used)
    end
end
%%
% Caculate intertia about weights own cg.
Jwdo = (1 / 2) * mwd * rwdo^2;

%%
% Combine drive inertia:
Jd = Jdd + mwd * rwd^2 + Jwdo;

%%
% Calculate Load inertia
if mwl < 0.81
    if mwl > 0.39
        rwlo = 0.016; % (smaller brass weight used)
    end
end

if mwl < 2.1
    if mwl > 0.9
        rwlo = 0.025; % (larger brass weight used)
    end
end

%%
% Calculate inertia about weights own cg.
Jwlo = (1 / 2) * mwl * rwlo^2;

%%
%Combine load inertia
Jl = Jdl + mwl * rwl^2 + Jwlo;

%%
% Build transfer function and state space models
if Config == 1  % Drive disk only:
    %Transfer function:
    N = khw;
    D = [Jd c1 0];
    %State Space model:
    A1 = [0 1];
    A2 = [0 -c1/Jd];
    Aol = [A1; A2];
    B = [0 khw/Jd]';
    C = [1 0]; %Theta1 output
end

if Config == 2  % Drive & load disk -- rigid drive train
    if Out == 1 %Encode %1 output
        Jr = Jd + Jp * grprime^(-2) + Jl * gr^(-2);  %Reflected inertia at drive
        cr = c1 + c2 * gr^(-2);  %Relected damping to drive
        %Transfer Function:
        N = khw;
        D = [Jr cr 0];
        %State space model:
        A1 = [0 1];
        A2 = [0 -cr/Jr];
        Aol = [A1; A2];
        B = [0 khw/Jr]';
        C = [1 0];
    end
    if Out == 2 %Encode %2 output
        Jr = Jd * gr^2 + Jp * (gr / grprime)^2 + Jl; %Reflected inertia at load
        cr = c1 * gr^2 + c2;     %Reflected dampling at load
        %Transfer Function
        N = khw * gr;
        D = [Jr cr 0];
        %State space model:
        A1 = [0 1];
        A2 = [0 -cr/Jr];
        Aol = [A1; A2];
        B = [0 khw*gr/Jr]';
        C = [1 0];
    end
end

if Config == 3  %Drive & load disks -- flexible drive train
    Jdstr = Jd + Jp * grprime^(-2); %Pulley inertias combined with drive
    
    %Transfer Function
    %The following do not include the coupled damping c12
    N1 = khw * [Jl c2 k];
    N2 = khw * k / gr;
    D = [Jdstr*Jl (c2*Jdstr+c1*Jl) (k*(Jdstr+Jl/gr^2)+c1*c2) (k*(c1+c2/gr^2)) 0];
    
    %The following include c12
    %N1 = khw*[Jl c2+c12 k];
    %N2 = khw/gr*[c12 k];
    %D = [Jdstr*Jl (c2*Jdstr+c1*Jl+(Jdstr+Jl/gr^2)*c12) (k*(Jdstr+Jl/gr^2)+c1*c2+(c1+c2/gr^2)*c12) (k*(c1+c2/gr^2)) 0];
    
    
    %State space model
    A1 = [0 1 0 0];
    A3 = [0 0 0 1];
    %The following do not include the coupled damping c12
    A2 = [-k/Jdstr/gr^2 -c1/Jdstr k/Jdstr/gr 0];
    A4 = [k/Jl/gr 0 -k/Jl -c2/Jl];
    
    %The following include c12
    %A2 = [-k/gr^2/Jdstr -(c1+c12/gr^2)/Jdstr k/gr/Jdstr c12/gr/Jdstr];   
    %A4 = [k/gr/Jl c12/gr/Jl -k/Jl -(c2+c12)/Jl];
    
    Aol = [A1; A2; A3; A4];
    B = [0 khw/Jdstr 0 0]';
    if Out == 1 %Encode %1 output
        C = [1 0 0 0];
    end
    if Out == 2 %Encode %2 output
        C = [0 0 1 0];
    end
end

% End of construct the plant 

%% 2. Construct transfer function
% Write transfer function directly
s = tf('s');

if Config == 1 || Config == 2
    num = N / D(1);
    den(2) = D(2) / D(1);
    den(3) = 0;
    den(1) = 1;
    H = tf(num, den);
end

if Config == 3
    num = N2 / D(1);
    den(2) = D(2) / D(1);
    den(3) = D(3) / D(1);
    den(4) = D(4) / D(1);
    den(1) = 1;
    den(5) = 0;
    H = tf(num, den);
end

% Write the system using A, B, C, and D
sys = ss(Aol, B, C, 0);



%% 3. Find Controllable, Observable, and Jordam Canonical Form
% Controllable and Observable form can be found by hand 
[Ac, Bc, Cc, Dc, P1] = canon(Aol, B, C, 0, 'companion'); % This is observable form

% Jorndan form
[V,J] = jordan(Aol);
[AJ, BJ, CJ, DJ, P2] = canon(Aol, B, C, 0);


%% 4. Find impulse response and step response
% Find close loop system

Hcl = 1/(1+sys); % or use command Hcl = feedback(sys, 1);


% % Impulse response
% 
% figure
% impulse(Hcl);
% grid;
% 
% % Step response
% 
% figure
% stepplot(Hcl);
% grid;

%% 5. Plot the Bode plot of the uncompensated system as well as root-locus of open loop system

% % Bode plot
% figure
% bode(sys);
% grid;
% 
% 
% % Root locus
% 
% figure
% rlocus(sys);
% grid;


%% 6. Design a lead-lag or a PID controller to meet certain design specifications (of your own choice). 
% Try to include both transient was well as steady state characteristics.
%
% Here is the design specification:
% $t_s$ <= 2s;
% P.O. < 5%; 
% $e_{ss}$ due to a step input is zero

%% 
%
% $$t_s \le 2 s $$  --> $\xi \omega \ge 2$ --> $\omega \ge \frac{2}{\xi} =2.86$; 
%
% $$P.O. \le 5 $$ --> $\theta < 45$ and $\xi = 0.707$; 
%
% $e_{ss}$ due to a step input is zero  --> put a 1/s in the forward path
%
% Design lead-lag compensator by using root-locus.
%
% The controller Gc is written as:

%%
% 
% $G_c = 0.012653 \frac{(1+0.56s)}{(1+0.23s)}$
%

if Config ==2
    Gc = 0.012134 * (1+0.56*s)/(1+0.25*s);
end

if Config == 3
    Gc = 0.012134 * (1+0.56*s)*(1+0.00091*s+(0.024*s)^2)/(1+0.25*s);
end

sys_controlled_op = Gc*H;
sys_controlled = Gc*sys/(1+Gc*sys);

%% 7. Obtain the step response, square wave and sinusoidal resposnes.

t = 0:0.01:10;

%% 
% Plot controlled step response

% figure
% step(sys_controlled);
% grid;

%%
% Square wave response;

[squareWave, tt] = gensig('square', 4, 10, 0.01);
% figure
% ysq = lsim(sys_controlled, squareWave, tt);
% plot(tt, squareWave, '--k', tt, ysq, '-r', 'LineWidth', 2);
% xlabel('Time (sec)');
% ylabel('Square wave response');
% legend('Input', 'Square wave response');
% grid;


%%
% Sinusoidal wave response;

u = 0.01*sin(t); % input

% figure
% ysq = lsim(sys_controlled, u, t);
% plot(t, u, '--k', t, ysq, '-r', 'LineWidth', 2);
% xlabel('Time (sec)');
% ylabel('Sinusoidal wave response');
% legend('Input = 0.01sin(t)', 'Sinusoidal response');
% grid;

%% 9. Design a full state feedback control to meet the design specification indicated in (6)
% $\omega \ge 2.86$ 
%  and $\xi = 0.707$.
% 
% So characteristic equation is $s^2 + 2\xi \omega s + \omega^2 = s^2 + 4.04404s + 8.1796$
% ** This is for second order system, so Config == 2

xi = 0.707;
omega = 2.86;


% Find pole

if Config == 2
    poles = [(-2*xi*omega+sqrt((2*xi*omega)^2-4*omega^2))/2 (-2*xi*omega-sqrt((2*xi*omega)^2-4*omega^2))/2];
end

if Config ==3
    % To design a controller 4th order system, firstly satisfy the 
    % requirements of 2nd order system, then move other pole far way from 
    % jw axis to make them less dominent.
    % So take two poles same as 2nd order system, then 
    poles = [-20 -15 (-2*xi*omega+sqrt((2*xi*omega)^2-4*omega^2))/2 (-2*xi*omega-sqrt((2*xi*omega)^2-4*omega^2))/2];
end

% Find K
    K = place(Aol, B, poles);
    Nbar = rscale(sys, K);
    
% Contructed controlled system by state feedback
    sys_cl = ss(Aol-B*K, B*Nbar, C, 0);

%% 10. Obtain the step, square wave and sinusoidal responses with arbitrary initial conditions and compare the results with those in 7.

%%
% Step response

% figure
% stepplot(sys_cl);
% grid;

%%
% Square wave response-- use lsim command

% figure
% ysq2 = lsim(sys_cl, squareWave, tt);
% plot(tt, squareWave, '--k', tt, ysq2, '-r', 'LineWidth', 2);
% xlabel('Time (sec)');
% ylabel('Square wave response');
% legend('Input', 'Square wave response');
% title('Square wave response using state feedback control');
% grid;
%%
% Sinusoidal wave response;
% 
% figure
% y1 = lsim(sys_cl, u, t);
% plot(t, u, '--k', t, y1, '-r', 'LineWidth', 2);
% xlabel('Time (sec)');
% ylabel('Sinusoidal wave response');
% legend('Input = 0.01sin(t)', 'Sinusoidal response');
% title('Sinusoidal wave response using state feedback control');
% grid;

%% 11. Design a full-order and a reduced-order observer and obtain step and sinusoidal response.

% --------------------
% Full-order observer -- with 5 times faster than designed controller

L = place(Aol', C', 5*poles);

At = [Aol-B*K B*K;
    zeros(size(Aol)) Aol-L'*C];
Bt = [B*Nbar; zeros(size(B))];
Ct = [C zeros(size(C))];

sys_est = ss(At, Bt, Ct, 0);

%% 
% Step response

% figure
% stepplot(sys_est); grid;

%%
% Sinusoidal wave response;

% figure
% y2 = lsim(sys_est, u, t);
% plot(t, u, '--k', t, y2, '-r', 'LineWidth', 2);
% xlabel('Time (sec)');
% ylabel('Sinusoidal wave response');
% legend('Input = 0.01sin(t)', 'Sinusoidal response');
% title('Sinusoidal wave response of full-order observer');
% grid;

% --------------------------
% Reduced-order observer

% Define requirements for reduced-order observer
% Define observer poles

if Config == 2
    % Reduce to 1st order system, with 3\tau < 2 --> \tau = 0.5 and \xi = 1
    % 5 times faster than controller.
    desiredPoles = -3*5;
    Aaa = [0];
    Aab = [1];
    Aba = [0];
    Abb = [-1.7820];
    
    L_redu = acker(Abb', Aab', desiredPoles')';
    A4reduObs = Abb - (L_redu * Aab);
    B4reduObs = eye(1);
    C4reduObs = eye(1);
end

if Config == 3
    % Reduce to 3rd order system, with desired poles of state feedback
    % controller except -20. 5 times faster than controller
    
    desiredPoles = 5*[-15 (-2*xi*omega+sqrt((2*xi*omega)^2-4*omega^2))/2 (-2*xi*omega-sqrt((2*xi*omega)^2-4*omega^2))/2];
    Aaa = Aol(1);
    Aab = Aol(1, 2:4);
    Aba = Aol(2:4, 1);
    Abb = Aol(2:4, 2:4);
    
    L_redu = acker(Abb', Aab', desiredPoles')';
    A4reduObs = Abb - (L_redu * Aab);
    B4reduObs = B(2:4);
    C4reduObs = C(2:4);
end

%Reduced-order estimator
sys_reduEst = ss(A4reduObs, B4reduObs, C4reduObs, 0); 

%% 
% Step response
% 
% figure
% stepplot(sys_reduEst); grid;

%%
% Sinusoidal wave response;
% 
% figure
% y3 = lsim(sys_reduEst, u, t);
% plot(t, u, '--k', t, y3, '-r', 'LineWidth', 2);
% xlabel('Time (sec)');
% ylabel('Sinusoidal wave response');
% legend('Input = 0.01sin(t)', 'Sinusoidal response');
% title('Sinusoidal wave response of reduced-order observer');
% grid;


%% 12. Find the transfer function of the oberserver and controller. What type of controller do you get?

Gec_D = (s*eye(length(Aol)) - Aol + B*K + L'*C);  % D(s) of Gec
%End of Plant.m 