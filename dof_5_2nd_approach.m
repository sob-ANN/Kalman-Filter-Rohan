% relative displacement of a mass considered w.r.t ground
add = 0;
T = 20 + add;
dt = 0.02;
t = 0:dt:T;
n = 5; % dof

% mass
m1 = 1;         % in order 1e+0
m2 = 1;
m3 = 1;
m4 = 1;
m5 = 1;

% stiffness
k1 = 10^2;       % in order 1e+2
k2 = 10^2;
k3 = 10^2;
k4 = 10^2;         % natural freq < sampling freq ( keep atleast 4 times )
k5 = 10^2;

% damping coefficient
c1 = 0.2;
c2 = 0.2;
c3 = 0.2;       % in order 1e-1
c4 = 0.2;
c5 = 0.2;

g_acc = load('ground_acc.mat','T').T';
g_acc = [g_acc zeros(1,add*50)];
u = zeros(5,size(t,2));
%u = 10*randn(5,size(t,2));
%g_acc = 0.0001*randn(1,length(t));
%u = [u; g_acc];
M = diag([m1,m2,m3,m4,m5]);
K = [k1 -k1 0 0 0;...
    -k1 k1+k2 -k2 0 0;...
    0 -k2 k2+k3 -k3 0;...
    0 0 -k3 k3+k4 -k4;...
    0 0 0 -k4 k4+k5];
C = [c1 -c1 0 0 0;
    -c1 c1+c2 -c2 0 0;...
    0 -c2 c2+c3 -c3 0;...
    0 0 -c3 c3+c4 -c4;...
    0 0 0 -c4 c4+c5];
% lets check natural frequency and zeta/damping ratio before proceeding
% further
[phi, L] = eig(K,M);
omg = sqrt(diag(L));
zeta = diag(phi'*C*phi)./omg/2;
max(zeta)
natural_freq = omg/(2*pi);
% Ac_ref = [zeros(n,n) eye(n);...
%     -k1/m1 k1/m1 0 0 0 -c1/m1 c1/m1 0 0 0;...
%     k1/m2 -(k1+k2)/m2 k2/m2 0 0 c1/m2 -(c1+c2)/m2 c2/m2 0 0 ;...
%     0 k2/m3 -(k2+k3)/m3 k3/m3 0 0 c2/m3 -(c2+c3)/m3 c3/m3 0 ;...
%     0 0 k3/m4 -(k3+k4)/m4 k4/m4 0 0 c3/m4 -(c3+c4)/m4 c4/m4 ;...
%     0 0 0 k4/m5 -(k4+k5)/m5 0 0 0 c4/m5 -(c4+c5)/m5];
Ac = [zeros(n,n) eye(n);...
    -M\K -M\C];
Bc = [zeros(n,n);...
    inv(M)];
Ec = [zeros(n,1);-ones(n,1)];
C = Ac(n+1:2*n,:);
D = Bc(n+1:2*n,:);
Ad = expm(Ac*dt);
Bd = Ac\(Ad - eye(2*n))*Bc;
Ed = Ac\(Ad - eye(2*n))*Ec;

% simulate
x0 = zeros(2*n,1);
x = x0;
x_store = zeros(2*n,size(t,2));
y_store = zeros(n,size(t,2));

for i = 1:length(t)
    x = Ad*x + Bd*u(:,i) + Ed*g_acc(:,i);
    y = C*x + D*u(:,i); 
    x_store(:,i) = x;
    y_store(:,i) = y;
end
plot(t,y_store(1,:))
ylabel('acceleration')
xlabel('time in seconds')