
function [INN,K,Pp,Pf,Yp,Yf,Xp,Xf] = Kalm(A,B,C,G,U,Y,Q,R)

% function [INN,K,Pp,Pf,Yp,Yf] = kalm(A,B,C,G,U,Y,Q,R)
%
% KALM implementa il filtro di Kalman per il sistema seguente :
%
%
%		x(k+1) = A x(k) + B u(k) + G w(k)
%		             
%		y(k) = C x(k) + v(k)
%		
% dove w(k) e v(k) sono rumori bianchi incorrelati con matrici di covarianza
% rispettivamente Q ed R. 
%    
% Uscite : INN -> innovazioni
%	  	   K   -> guadagno del filtro di Kalman a regime
%          Pp   -> matrice di covarianza dell'errore di stima a regime (predizione)
%          Pf   -> matrice di covarianza dell'errore di stima a regime (filtraggio)
%          Yp  -> uscite predette
%          Yf  -> uscite filtrate
                    
[L,ncU] = size(U);
[LY,ncY] = size(Y);
[n,ncA] = size(A);
[nrB,r] = size(B);
[m,ncC] = size(C);
[nrQ,ncQ] = size(Q);
[nrR,ncR] = size(R);

% Test sulla congruenza dei dati di ingresso

if L ~= LY,  error('Dimensional incongruence -> U,Y'); end
if n ~= nrB, error('Dimensional incongruence -> A,B'); end
if n ~= ncC, error('Dimensional incongruence -> A,C'); end
if ncU ~= r, error('Dimensional incongruence -> U,B'); end
if ncY ~= m, error('Dimensional incongruence -> Y,C'); end
%if ncQ ~= r, error('Dimensional incongruence -> U,Q'); end
if ncR ~= m, error('Dimensional incongruence -> Y,R'); end

% Initial state and covariance matrix of the estimate

INN = zeros(L,m);
Yp = zeros(L,m);
Yf = zeros(L,m);
Xf = zeros(L,n);
%Xkp = pinv([C;C*A])*[Y(1,:)';Y(2,:)'];
Xkp=Y(1:n);
%Xkp=zeros(n,1);
%Xp(1,:)=Xkp';
Pkp = eye(n);
Yp(1,:) = (C*Xkp)';
INN(1,:) = Y(1,:) - Yp(1,:);
Kk = Pkp*C'*inv(C*Pkp*C' + R);
Xkf = Xkp + Kk*INN(1,:)';
Yf(1,:) = (C*Xkf)';
Xf(1,:) = Xkf';
Pkf=Pkp - Kk*C*Pkp;
K=Kk';

for k = 1:L-1
  Xkp = A*Xkf + B*U(k,:)';
  Xp(k+1,:)=Xkp';
  Pkp = A*Pkf*A' + G*Q*G';
  Yp(k+1,:) = (C*Xkp)';
  INN(k+1,:) = Y(k+1,:) - Yp(k+1,:);
  Kk = Pkp*C'*inv(C*Pkp*C' + R);
  K=[K;Kk'];
  Xkf = Xkp + Kk*INN(k+1,:)';
  Yf(k+1,:) = (C*Xkf)';  
  Xf(k+1,:) = Xkf';
  Pkf=Pkp - Kk*C*Pkp; 
end

Pp=Pkp;
Pf=Pkf;
%K=Kk;