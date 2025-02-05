function [p] = peak_LMS(x,L,debug)

if nargin < 3
    debug = 0;
end
alpha = 1;
x = x(:); %always column vector
N = size(x,1);

%L = ceil(N/2)-1; %maximum
%L = 200;   %if peaks are more frequent a lower L suffices
M = rand(L,N)+alpha;
for k = 1:L
    for i=k+2:N-k+1
        if (x(i-1)>x(i-k-1)) && (x(i-1)>x(i+k-1))
            M(k,i) = 0;
        end
    end
end

g          = sum(M,2);
[mm, lam]  = min(g);
Mi         = M(1:lam,:);
si         = std(Mi);
p          = find(si==0);


if debug
    figure(2);
    subplot(4,1,1);
    imagesc(M);
    subplot(4,1,2);
    plot(g);
    subplot(4,1,3);
    plot(si);
    subplot(4,1,4);
    plot(1:N,x,'b-');hold on;
    plot(p, x(p),'ro');hold off;
    pause
end
end

