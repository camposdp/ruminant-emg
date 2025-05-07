
close all
clc
clear



x = [0:.01:9.99];
y1 = normpdf(x,2,0.4);
y2 = normpdf(x,4,0.4);

figure
subplot(3,1,1)
hold
plot(x,y1)
plot(x,y2)
line([3 3],[0 1],'Color','k','LineStyle','--')
ylim([0 1])

ylabel('PDF')
xlabel('Característica')
legend('Classe 1','Classe 2','Discriminante')


y1 = normpdf(x,4,0.4);
y2 = normpdf(x,6,0.4);
subplot(3,1,2)
hold
plot(x,y1)
plot(x,y2)
line([3 3],[0 1],'Color','k','LineStyle','--')
ylim([0 1])

ylabel('PDF')
xlabel('Característica')
legend('Classe 1','Classe 2','Discriminante')

y1 = normpdf(x,4,1);
y2 = normpdf(x,6,1);
subplot(3,1,3)
hold
plot(x,y1)
plot(x,y2)
line([3 3],[0 1],'Color','k','LineStyle','--')
ylim([0 1])

ylabel('PDF')
xlabel('Característica')
legend('Classe 1','Classe 2','Discriminante')

x = [0:.01:99.99];

for i = 1:200
y1 = normpdf(x,2+i*(2/100),0.4);
y2 = normpdf(x,4+i*(2/100),0.4);

eu(i)= sum(y1(300:end))+sum(y2(1:300));

y1 = normpdf(x,2+i*(2/100),0.4+i*(0.6/100));
y2 = normpdf(x,4+i*(2/100),0.4+i*(0.6/100));

et(i) =  sum(y1(300:end))+sum(y2(1:300));

end


figure
subplot(1,2,1)
plot(eu)
subplot(1,2,2)

plot(et)
ylim([0 100])