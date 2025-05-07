


myColorMap = jet(256);
myColorMap(1,:) = 1;
x=zeros(25,10);

for i = 1:10
x(randi(25,5,1),i)=randi(100,5,1)./100;
end

x = x./max(abs(x));

pcolor(x)
colormap(myColorMap)