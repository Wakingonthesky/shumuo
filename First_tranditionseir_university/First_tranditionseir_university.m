clc
clear
close all;

S = 50000;
E = 0;
I = 1;
A = 0;
R = 0;
D = 0;


%疫苗接种率87%
%易感者变为感染者的速率0.33
%治愈率γ=0.154（每天每一百个中治愈15.4人）
%潜伏者转化为无症状感染者的比率alpha=1/4
%潜伏者转轻症为1/7
%潜伏者自愈0.2
%平均每人每天可以接触大概15人左右，即r=15
%不加任何干预下的死亡率为K1=0.0675
%无症状转轻症H=0.86
%病亡率K1=5.6%
%居家隔离占比80%
%老年人比例18.7%，患病率是成年人的5.4倍，死亡率3/1


N = S + I;

a = 1/4;  %潜伏者转化为感染者的速率
a1 = 0.2; %无症状感染者的治愈的概率
b = 1/7;  %潜伏者转轻症
cba = 0.8; %被接触感染的概率
d = 0.2;  %潜伏者自愈
mu = 1/7; %潜伏者转轻症
g = 0.154; %治愈率
k = 0.056; %死亡率
k1 = 0.0675;%不加任何干预下的死亡率
l = 0.99; %无症状感染者的治愈的速率
H = 0.86; %无症状转轻症
o = 0.187;%老年人比例

tspan = [0 150];
y0 = [S E A I R D];
[t, y] = ode45(@(t,y)odefun(t,y,a,a1,b,cba,H,g,k,N), tspan, y0);

plot(t,y(:,1),'b',t,y(:,2),'m',t,y(:,3),'r',t,y(:,4),'k',t,y(:,5),'c',t,y(:,6),'g')
xlabel('day')
ylabel('person')
legend('S','E','A','I','R','D')
title('First changed SEIR model(university)')

function dydt = odefun(t,y,a,a1,b,cba,H,g,k,N)
dydt = zeros(6,1);
dydt(1) = -cba*y(1)*y(4)/N-cba*y(1)*y(3)/N;
dydt(2) = cba*y(1)*y(4)/N+cba*y(1)*y(3)/N-b*y(2)-a*y(2);
dydt(3) = a*y(2)-H*y(3)-a1*y(3);
dydt(4) = b*y(2)+H*y(3)-g*y(4)-k*y(4);
dydt(5) = g*y(4)+a1*y(3);
dydt(6) = k*y(4)
end
