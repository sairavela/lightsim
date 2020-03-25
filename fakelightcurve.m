
x = randperm(100);lx = length(x);
y = zeros(1,1000);ly = length(y);
k  =200;
for i = 1:k:ly
    toss = rand>0.1;
    y(i) = -toss*0.02*(1+rand);
end
ly = length(y);


ker = [0:20 20:-1:0]/10;
x = randperm(ly);
off = zeros(size(x));
for N = 2:5
off = off + 1/N^2*sin(2*pi*N/ly*(1:ly)+rand*2*pi);
end
off = off/(200*rand+1);
z = off+0.985+(rand(size(y))/10+1).*conv(y,ker,'same')+randn(1,ly)/(100+50*rand);
z(x(4:10:end))=nan;
plot(z,'.')