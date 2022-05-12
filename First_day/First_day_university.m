clc
clear
close all;

S = 6000000;
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
%老年人比例18.7%，患病率是成年人的2倍，死亡率3/1


N = S + I;

a = 2/3;  %潜伏者转化为感染者的速率
a1 = 0.01;
b = 1/7;  %潜伏者转轻症
cba = 0.534; %被接触感染的概率概率
d = 0.2;  %潜伏者自愈
mu = 1/7; %潜伏者转轻症
g2 = 0.34;
 %死亡率
k1 = 0.225;
l = 0.99; %无症状感染者的治愈的
l1 = 1/5;
l2 = 1/3;
H = 0.86; %无症状转感染者
o = 0.187;%老年人比例
x1 = 1/3;
x2 = 2/3;
%x3 = 1/3;
v = 1/7;


y0 = [S E E1 A I I1 I2 R D];
tspan2 = [0 500];
[t, y1] = ode45(@(t,y)odefun(t,y,a,a1,b,cba,d,H,g2,k1,l1,12,x1,x2,N,v), tspan2, y0);
disp(sum(y1(:,4))+sum(y1(:,7)))
%%
tspan2 = [0 65.6];
[t, y2] = ode45(@(t,y)odefun(t,y,a,a1,b,cba,d,H,g2,k1,l1,12,x1,x2,N,v), tspan2, y0);

y1 = y2(end,:);

tspan2 = [65.6 500];
[t, y2] = ode45(@(t,y)odefun(t,y,a,a1,b,cba,d,H,g2,k1,l1,12,x1,x2,N,v), tspan2, y1);
disp(sum(y2(:,4))+sum(y2(:,7)))
%%
tspan3 = [0 78];
[t, y3] = ode45(@(t,y)odefun(t,y,a,a1,b,cba,d,H,g2,k1,l1,12,x1,x2,N,v), tspan3, y0);

y1 = y3(end,:);

tspan2 = [78 500];
[t, y3] = ode45(@(t,y)odefun(t,y,a,a1,b,cba,d,H,g2,k1,l1,12,x1,x2,N,v), tspan2, y1);
disp(sum(y3(:,4))+sum(y3(:,7)))
%%
tspan4 = [0 168];
[t, y4] = ode45(@(t,y)odefun(t,y,a,a1,b,cba,d,H,g2,k1,l1,12,x1,x2,N,v), tspan4, y0);

y1 = y4(end,:);

tspan2 = [168 500];
[t, y4] = ode45(@(t,y)odefun(t,y,a,a1,b,cba,d,H,g2,k1,l1,12,x1,x2,N,v), tspan2, y1);

disp(sum(y4(:,4))+sum(y4(:,7)))
%%
tspan5 = [200 500];
[t, y5] = ode45(@(t,y)odefun(t,y,a,a1,b,cba,d,H,g2,k1,l1,12,x1,x2,N,v), tspan5, y0);

y1 = y5(end,:);

tspan2 = [65.6 500];
[t, y5] = ode45(@(t,y)odefun(t,y,a,a1,b,cba,d,H,g2,k1,l1,12,x1,x2,N,v), tspan2, y1);
disp(sum(y5(:,4))+sum(y5(:,7)))
%%
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