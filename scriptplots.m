close all;
z_val=zeros(1,300);
y_val=zeros(1,300);
x_val=zeros(1,300);
w_val=zeros(1,300);
V_BS=0;
L=1*10^(-7);
W=3*10^(-7);
TEMP=27+273.15;

%%Subumbral en funcion de Vgs, Vds=0.5V
V_DS=0.2;
V_GS=linspace(-3,1,300);
for i=1:300
	[~,y_val(i),~,~]=nmos(V_DS,V_BS,V_GS(i),L,W,TEMP);
end
x_val=V_GS;
plot(x_val,y_val);
figure;

%%Subumbral en funcion de Vds, barrido de Vgs
V_DS=linspace(0,3,300);
V_GS=linspace(0,1.2,10);
hold on;
for j=1:10
    for i=1:300
        [z_val(i),y_val(i),~,w_val(i)]=nmos(V_DS(i),V_BS,V_GS(j),L,W,TEMP);
    end
x_val=V_DS;
plot(x_val,y_val);
end
hold off;
figure;
plot(x_val,z_val);
figure;
plot(x_val,w_val);

