%%Chaining Data under a Gaussian Distribution
% chaining invariant statistics 
% example code, sai ravela, (C) Copyright MIT, All rights reserved.

clear all;

N = 100; % Half length of buffer >16

H0 = zeros(2*N);
Hl = H0; Hl(1:N,1:N) = eye(N);Hl = Hl(1:N,:); %Left Half Operator
Hr = H0; Hr(N+1:end,N+1:end) = eye(N); Hr = Hr(N+1:end,:); %Right Half Operator

s = 20; % Correlation scale <=N/5
x = -N:N-1; % zero centered indices
y = 1/sqrt(2*pi)/s*exp(-(x+N).*(x+N)/2/(s.*s)); % Kernel
Cyy = toeplitz(y,y); % Covariance Function

[U,S] = eig(Cyy); % It is better if this is sorted for truncated modes -- Cyy is SPD.
S = real(sqrt(S)); % SQUARE roots
Cyy = U*(S*S')*U'; % ensure SPD
z =  U*S*randn(2*N,1); %Generate a random vector

for i = 1:200, % Iterations each iteration will add an extra N elements to z.
    lN = length(z);
    zthis = z(lN-2*N+1:lN); % The last 2N vector.
    z1 = U*S*randn(2*N,1); % Generate a random 2N vector
    
    %See ML class notes on "One Step at A Time"
    z2 = z1+Cyy*Hl'*pinv(Hl*Cyy*Hl')*(Hr*zthis-Hl*z1); %Demand that right half z1 and left half z2 match.
    
    subplot(121);   
    plot(1:lN,z);
    hold on
        plot(lN-N+1:lN+N,[z2],'g');
        plot(lN-N+1:lN,Hr*zthis,'r');
    hold off;%axis([1 lN+N -3 3]);
    
    z = [z; Hr*z2];
    
    subplot(122);
    plot(lN-2*N+1:lN,zthis);
    hold on;
        plot(lN-N+1:lN+N,[z2],'g');
        plot(lN-N+1:lN,Hr*zthis,'r');
        hold off;
        %axis([lN-2*N+1 lN+N -3 3]);
    pause(0.1);
end
