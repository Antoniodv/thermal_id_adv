% Linear Kalman filter
%
% Adaptation of file from:
% Version 1.0, May 16,2011
% Version 2.0, May 10,2013
% Written by Chong You
% https://sites.google.com/site/chongyou1987/
% chong.you1987@gmail.com
%
% Syntax:
%     [Xo,Po] = Extended_KF(f,g,Q,R,Z,Xi,Pi)
%
% State Equation:
%     X(n+1) = A * X(n) + B * U(n) + w(n)
%     where the state X has the dimension N-by-1
% Observation Equation:
%     Z(n) = C * X(n) + v(n)
%     where the observation y has the dimension M-by-1
%     w ~ N(0,Q) is gaussian noise with covariance Q
%     v ~ N(0,R) is gaussian noise with covariance R     
% Input:
%     A: matrix of state transition
%     B: matrix of control input
%     C: matrix of output update
%     Q: process noise covariance matrix, N-by-N
%     R: measurement noise covariance matrix, M-by-M -> assumed diagonal
%     with standard deviation meas_sigma
%     
%     Xi: "a priori" state estimate, N-by-1
%     Pi: "a priori" estimated state covariance, N-by-N
%     Uk: current input
%     zk: current measurements
% Output:
%     Xo: "a posteriori" state estimate, N-by-1
%     Po: "a posteriori" estimated state covariance, N-by-N
%
% Algorithm for Extended Kalman Filter:
% Linearize input functions f and g to get fy(state transition matrix)
% and H(observation matrix) for an ordinary Kalman Filter:
% State Equation:
%     X(n+1) = fy * X(n) + w(n)
% Observation Equation:
%     Z(n) = H * X(n) + v(n)
%
% 1. Xp = A * Xi + B * Ui           : One step projection, also provides 
%                                     linearization point
% 
% 2. 
%          d f    |
% fy = -----------|                 : Linearize state equation, fy is the
%          d X    |X=Xp               Jacobian of the process model
%       
% 
% 3.
%          d g    |
% H  = -----------|                 : Linearize observation equation, H is
%          d X    |X=Xp               the Jacobian of the measurement model
%             
%       
% 4. Pp = fy * Pi * fy' + Q         : Covariance of Xp
% 
% 5. K = Pp * H' * inv(H * Pp * H' + R): Kalman Gain
% 
% 6. Xo = Xp + K * (Z - g(Xp))      : Output state
% 
% 7. Po = [I - K * H] * Pp          : Covariance of Xo
	                                                                            
function [Xo,Po,zf,yi] = KalmanFilter(A, B, C, Q, R, Xi, Pi, Uk, zk)

% Predict
Xp = A * Xi + B * Uk;
Pp = A * Pi * A.' + Q;

% Update
yi = zk - C * Xp;
Sk = C * Pp * C.' + R;
K = Pp * C' / Sk;  % in reality it is: Pp * C' * inv(Sk);
Xo = Xp + K * yi;
zf = C * Xo;
Po = Pp - K * C * Pp;
yo = zk - C * Xo;
    
 