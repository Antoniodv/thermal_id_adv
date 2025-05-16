function [theta,sigv,sigw,res,condn,svd_ident] = misoarxbls(U,Y,n,num_add_eq,theta_init)

% Bilinear least squares estimate of MISO ARX models with additive output noise
% U = [U1 U2 ... UnU] input sequences
% Y output sequence
% n model order
% theta = [a1 ... an b1 ... bn]' estimated parameter vector
% sigv estimated variance of the additive noise

if (nargin < 4)
	num_add_eq = 10;
end
if (nargin < 5)
    theta_init = [];
end


[L,nU]=size(U);
% for i=1:nU
%     U(:,i)=U(:,i)-mean(U(:,i));
% end
L=length(Y);
% Y=Y-mean(Y);

N=L-n;

% Hankel matrices of input and output signals

Hy=hank(Y,n+1,L-n);
Hu=[];
for i=1:nU
    Hui=hank(U(:,i),n+1,L-n);
    Hui(:,1)=[];
    Hu=[Hu Hui];
end

H=[-Hy Hu];

% Computation of R and r
yy=-H(:,1);
Hls=H;Hls(:,1)=[];
R=Hls'*Hls/N;
r=Hls'*yy/N;

% Augmented Hankel matrices of input signals

q=num_add_eq;  % number of additional equations (d=n+q)

Zy=hank(Y,n+q+1,L-n-q);
Zu=[];
HH=-Zy(:,1:n+1);
for i=1:nU
    Zui=hank(U(:,i),n+q+1,L-n-q);
    Zui(:,1)=[];
    Zu=[Zu Zui];
    HH=[HH Zui(:,1:n)];
end

SIGMAh=Zu'*HH/(L-n-q);

% Computation of Rd and rd

Rd=SIGMAh;
Rd(:,1)=[]; 
rd=SIGMAh(:,1);

rho=[r;rd];
Rbar=[R;-Rd];
[n1,n2]=size(R);
[nh1,nh2]=size(Rd);
D=zeros(n1,n2);
D(1:n,1:n)=eye(n);
J=[D;zeros(nh1,nh2)];

thetaLS=pinv(R)*r; % least squares estimate
thetak=zeros(length(thetaLS),1); % 
if isempty(theta_init)
    thetakp1=thetaLS; % algorithm inizialization
else
    thetakp1=theta_init; % algorithm inizialization
end

epsilon = 1e-3;  % convergence threshold

k=0;
while norm(thetakp1-thetak)/norm(thetakp1) > epsilon 
thetak = thetakp1;
ak=thetak(1:n);
sigvk=(thetak'*J'*(Rbar*thetak-rho))/(ak'*ak);
Jk=J;Jk(1:n,1:n)=eye(n)*sigvk;
thetakp1=pinv(Rbar-Jk)*rho;
k=k+1;
end

theta=thetakp1;
sigv=sigvk;
res = H*[1;theta];

% Calculate the driving noise variance
sigmaT = mean(Y.^2);
sigw = sigmaT - r' * theta - sigv;

% Output the conditioning number
svd_ident = svd(Rbar-Jk);
condn = cond(Rbar-Jk);

end

