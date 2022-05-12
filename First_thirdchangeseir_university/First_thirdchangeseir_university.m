clc
clear
close all;

S = 50000;
E = 0;
E1 = 0;
I = 1;
I1 = 0;
I2 = 0;
I3 = 0;
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

a = 2/3;  %潜伏者转化为感染者的速率
a1 = 0.01;
b = 1/7;  %潜伏者转轻症
cba = 0.8; %被接触感染的概率概率
d = 0.2;  %潜伏者自愈
mu = 1/7; %潜伏者转轻症
%g = 0.154; %治愈率
g1 = 0.154;
g2 = 0.34;
%g3 = 0.154;
k = 0.5; %死亡率
k1 = 0.0057;
%k3 = 0.056;
l = 0.99; %无症状感染者的治愈的
l1 = 1/5;
l2 = 1/3;
H = 0.86; %无症状转感染者
o = 0.187;%老年人比例
x1 = 1/3;
x2 = 2/3;
%x3 = 1/3;
v = 1/7;

tspan = [0 40];
y0 = [S E E1 A I I1 I2 R D];
[t, y] = ode45(@(t,y)odefun(t,y,a,a1,b,cba,d,H,g2,k1,l1,12,x1,x2,N,v), tspan, y0);
hold on;

ha = 2/3;  %潜伏者转化为感染者的速率
ha1 = 0.01;
hb = 1/7;  %潜伏者转轻症
hcba = 0.5; %被接触感染的概率概率
hd = 0.2;  %潜伏者自愈
hmu = 1/7; %潜伏者转轻症
%g = 0.154; %治愈率
hg1 = 0.154;
hg2 = 0.34;
%g3 = 0.154;
hk = 0.5; %死亡率
hk1 = 0.0057;
%k3 = 0.056;
hl = 0.99; %无症状感染者的治愈的
hl1 = 1/5;
hl2 = 1/3;
hH = 0.86; %无症状转感染者
ho = 0.187;%老年人比例
hx1 = 1/3;
hx2 = 2/3;
%x3 = 1/3;
hv = 1/7;

HS = y(end,1);
HE = y(end,2);
HE1 = y(end,3);
HI = y(end,4);
HI1 = y(end,5);
HI2 = y(end,6);
HA = y(end,7);
HR = y(end,8);
HD = y(end,9);


tspan = [40 150];
y0 = [HS HE HE1 HA HI HI1 HI2 HR HD];
[t1, y1] = ode45(@(t1,y1)odefun1(t1,y1,a,a1,b,cba,d,H,g2,k1,l1,12,x1,x2,N,v), tspan, y0);

plot(t1,y1(:,1),'b',t1,y1(:,2),'m',t1,y1(:,4),'r',t1,y1(:,6),'g',t1,y1(:,7),'k',t1,y1(:,8),'c',t1,y1(:,9))
xlabel('day')
ylabel('person')
legend('S','E','A','I1','I2','R','D')
title('First changed SEIR model(university)')


function dydt = odefun(t,y,a,a1,b,cba,d,H,g2,k1,l1,l2,x1,x2,N,v)
dydt = zeros(9,1);
dydt(1) = -cba*y(1)*(y(6)+y(7))/N-cba*y(1)*y(4)/N+v*y(8)+v*y(3); %s
dydt(2) = cba*y(1)*(y(6)+y(7))/N+cba*y(1)*y(4)/N-b*y(2)-a*y(2)-d*y(2); %e
dydt(3) = d*y(2)-v*y(3); %e1
dydt(4) = a*y(2)-H*y(4)-a1*y(4); %a
dydt(5) = b*y(2)+H*y(4)-(x1+x2)*y(5); %ri
dydt(6) = x1*y(5)+l1*y(7)-l2*y(6)-k1*y(6); %i1
dydt(7) = x2*y(5)-l1*y(7)+l2*y(6)-g2*y(7); %i2
%dydt(8) = x3*y(5)-g3*y(8)-k3*y(8); %i3
dydt(8) = g2*y(7)+a1*y(4)-v*y(8); %r
dydt(9) = k1*y(6); %d
end

function dydt = odefun1(t1,y1,a,a1,b,cba,d,H,g2,k1,l1,l2,x1,x2,N,v)
dydt = zeros(9,1);
dydt(1) = -cba*y1(1)*(y1(6)+y1(7))/N-cba*y1(1)*y1(4)/N+v*y1(8)+v*y1(3); %s
dydt(2) = cba*y1(1)*(y1(6)+y1(7))/N+cba*y1(1)*y1(4)/N-b*y1(2)-a*y1(2)-d*y1(2); %e
dydt(3) = d*y1(2)-v*y1(3); %e1
dydt(4) = a*y1(2)-H*y1(4)-a1*y1(4); %a
dydt(5) = b*y1(2)+H*y1(4)-(x1+x2)*y1(5); %ri
dydt(6) = x1*y1(5)+l1*y1(7)-l2*y1(6)-k1*y1(6); %i1
dydt(7) = x2*y1(5)-l1*y1(7)+l2*y1(6)-g2*y1(7); %i2
%dydt(8) = x3*y(5)-g3*y(8)-k3*y(8); %i3
dydt(8) = g2*y1(7)+a1*y1(4)-v*y1(8); %r
dydt(9) = k1*y1(6); %d
end