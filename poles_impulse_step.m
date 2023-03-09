function [y_impulse,t_impulse, y_step, t_step] = poles_impulse_step(p12_real,p12_im,p3)
% clear all; 
% close all;

syms s a b l g Kp Ki    % define symbolic variables

Hvtheta = -s/l/(s^2-g/l);       % TF from velocity to angle of pendulum

K = Kp + Ki/s;                  % TF of the PI angle controller
M = a*b/(s+a);                   % TF of motor (1st order model) 
%  M = 1;                       % TF without motor
%
%  
% closed loop transfer function from disturbance d(t)totheta(t)
Hcloop = 1/(1-Hvtheta*M*K);    % use this for no motor feedback

pretty(simplify(Hcloop));       % to display the total transfer function

% Substitute parameters and solve
% system parameters
g = 9.81;
l = 0.436;   %effective length 
a = 13.83;           %nominal motor parameters
b = 0.0036;        %nominal motor parameters

Hcloop_sub = subs(Hcloop); % sub parameter values into Hcloop

% p12_real = -3.9
% p12_im = -2*pi
% wn = sqrt(p12_real^2 + p12_im^2)


% p12_im = -1*pi
% wn = 4.8
% p12_real = sqrt(wn^2 - p12_im^2)


% specify locations of the target poles,
% choose # based on order of Htot denominator
% e.g., want some oscillations, want fast decay, etc. 
p1 = p12_real + p12_im*i;    % dominant pole pair
p2 = p12_real - p12_im*i;   % dominant pole pair 
% p3 = -1


% target characteristic polynomial
% if motor model (TF) is added, order of polynomial will increase
tgt_char_poly = (s-p1)*(s-p2)*(s-p3);
npoly = 3;


% get the denominator from Hcloop_sub
[n d] = numden(Hcloop_sub);

% find the coefficients of the denominator polynomial TF
coeffs_denom = coeffs(d, s);

% divide though the coefficient of the highest power term
coeffs_denom = coeffs(d, s)/(coeffs_denom(end));

% find coefficients of the target charecteristic polynomial
coeffs_tgt = coeffs(tgt_char_poly, s);

% solve the system of equations setting the coefficients of the
% polynomial in the target to the actual polynomials
solutions = solve(coeffs_denom(1:npoly-1) == coeffs_tgt(1:npoly-1),  Kp, Ki);

% display the solutions as double precision numbers
Kp = double(solutions.Kp)
Ki = double(solutions.Ki)

% reorder coefficients for the check polynomial 
for ii = 1:length(coeffs_denom)
    chk_coeffs_denom(ii) = coeffs_denom(length(coeffs_denom) + 1 - ii);
end
closed_loop_poles = vpa (roots(subs(chk_coeffs_denom)), npoly );


% Plot impulse response of closed-loop system
    TFstring = char(subs(Hcloop));
    % Define 's' as transfer function variable
    s = tf('s');
    % Evaluate the expression
    eval(['TFH = ',TFstring]);
%     figure (1)
    [y_impulse,t_impulse] = impulse(TFH);   %plot the impulse reponse
%     hold on
%     figure(2)
    [y_step,t_step] = step(TFH);       %plot the step response
end

