function [V_thres,I_drain,I_leak,V_early] = nmos(V_DS,V_BS,V_GS,L,W,TEMP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 Constantes                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_b=1.3806488.*10.^(-23); %--Constante de Boltzmann.
q=1.60217657.*10.^(-19); %--Carga elemental de un electrón.
n_ch=2.49865.*10.^(17); %--Concentración de dopado del substrato.
n_d=1.*10.^(20); %--Concentración de dopado de difusiones.
n_i=1.45.*10.^(10); %--Concentración de portadores intrínsecos.
e_si=11.68; %--Permitividad relativa del silicio.
e_ox=3.9; %--Permitividad relativa del óxido. (SiO2)
e_0=8.8541878176.*10.^(-12); %--Permitividad del vacío.
t_ox=2.80.*10.^(-9); %--Espesor del óxido.
u_0=542.8809244; %--Movilidad de electrón en vacío.
E_0=0.67*10.^(6); %--Campo eléctrico experimentado en estado estable.
v=1.6; %--Parámetro empírico de corrección para cálculo de movilidad efectiva.
velsat=113760; %--Velocidad máxima de portadores de carga en canal.
v_th0=0.121512; %--Tensión de umbral nominal para transistor de canal largo sin efecto de cuerpo.
xj=6.35926.*10.^(-08); %--Profundidad de difusiones.
k_1=0.28925; %--Parámetro empírico de ajuste de efecto de cuerpo para cálculo de efecto de carga de substrato (Bulk charge effect).
a_0=1.99651; %--Parámetro empírico de ajuste de longitud de canal para cálculo de efecto de carga de substrato.
a_gs=1.78401; %--Parámetro empírico de ajuste de tensión de compuerta para cálculo de efecto de carga de substrato.
b_0=1.33242.*10.^(-06); %--Parámetro empírico de ajuste de canal angosto para cálculo de efecto de carga de substrato.
b_1=9.43998.*10.^(-07); %--Parámetro empírico de offset de ancho de canal para cálculo de efecto de carga de substrato.
keta=-0.0770244; %--Parámetro empírico de ajuste de tensión de substrato para cálculo de efecto de carga de substrato.
v_off=-0.0833823; %--Offset empírico de tensión efectiva en región de subumbral.
n=2.32861; %--Pendiente de subumbral.
k_i=1; %--Parámetro empírico de ajuste de corriente de fuga entre drenador y substrato.
v_i=30; %--Parámetro empírico de ajuste de corriente de fuga entre drenador y substrato.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 Expresiones                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_t=k_b*TEMP/q; %--Voltaje térmico.
v_bi=v_t*log(n_ch*n_d/(n_i).^2); %--Potencial intrínseco en la frontera de la zona de agotamiento.
phi_s=2*v_t*log(n_ch/n_i); %--Potencial de superficie del substrato.
x_dep=sqrt(2*e_si*e_0*(phi_s-V_BS)/(q*n_ch)); %--Profundidad promedio de zona de agotamiento.
lambda=sqrt(e_si*t_ox*x_dep/(e_ox)); %--Longitud característica del canal.
lambda=1.5*10.^(-8);
c_ox=e_ox*e_0/t_ox*10.^(-4); %--Capacitancia del óxido por unidad de área.
gamma=sqrt(2*q*e_si*n_ch)/c_ox; %Coeficiente de efecto de cuerpo de canal largo.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%              Tensión de umbral                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_th_long=v_th0+gamma*(sqrt(phi_s-V_BS)-sqrt(phi_s)); %--Tensión de umbral de canal largo.
theta_vth=exp(-L/(2*lambda))+2*exp(-L/lambda); %--Dependencia de tensión de umbral con largo de canal.
delta_vth=theta_vth*(2*(v_bi-phi_s)+V_DS); %--Reducción de Vth debido a efecto de canal corto.
narrow_vth=3*pi*(t_ox/W)*phi_s; %--Aumento de Vth debido a efecto de canal angosto.
v_th=v_th_long-delta_vth+narrow_vth; %--Tensión de umbral de canal corto.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%          Movilidad y campo eléctrico           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_eff=(V_GS-v_th)/(6*t_ox); %--Campo eléctrico promedio experimentado por los portadores en el canal.
u_eff=u_0/(1+(E_eff/E_0).^v); %--Movilidad efectiva de los electrones en el canal.
E_sat=2*velsat/(u_0*10.^(-4)); %--Campo eléctrico máximo experimentado por los portadores en el canal.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Coeficiente de efecto de carga de substrato  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_bulk=(1+k_1/(2*sqrt(phi_s-V_BS))*(a_0*L/(L+2*sqrt(xj*x_dep))*(1-a_gs*(V_GS-v_th)*(L/(L+2*sqrt(xj*x_dep))).^2)+(b_0/(W+b_1))))*(1/(1+keta*V_BS)); %--Coeficiente de efecto de carga de substrato
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%             Corriente de drenador              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_ds_lin=u_0*c_ox*(W/L)*(1/(1+(V_DS/(E_sat*L))))*((V_GS-v_th)-a_bulk*V_DS/2)*V_DS; %--Corriente de drenador en región lineal.
v_sat=(E_sat*L*(V_GS-v_th))/(a_bulk*E_sat*L+(V_GS-v_th)); %--Tensión de saturación del transistor.
i_sat0=u_0*c_ox*(W/L)*(1/(1+(v_sat/(E_sat*L))))*((V_GS-v_th)-a_bulk*v_sat/2)*v_sat; %--Corriente de saturación sin resistencia de salida.
v_a=(a_bulk*E_sat*L+(V_GS-v_th))*(V_DS-v_sat)/(a_bulk*E_sat*lambda); %--Voltaje de Early, únicamente se considera CLM.
i_ds_sat=i_sat0*(1+(V_DS-v_sat)/10.5); %--Corriente de saturación con resistencia de salida.
i_sub0=u_0*(W/L)*sqrt((q*e_si*e_0.*10.^(-2)*n_ch)/(2*phi_s))*v_t.^2; %--Corriente pre-exponencial de subumbral.
i_ds_sub=i_sub0*(1-exp(-V_DS/v_t))*exp((V_GS-v_th-v_off)/(n*v_t)); %--Corriente de subumbral del transistor.
if V_GS < v_th
		i_ds=i_ds_sub;
elseif V_GS >= v_th
	if V_DS >= v_sat
		i_ds=i_ds_sat;
	else
		i_ds=i_ds_lin;
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%               Corriente de fuga                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_db=i_ds*k_i*(V_DS-v_sat)*exp(-(v_i)/(V_DS-v_sat)); %--Corriente de fuga entre drenador y substrato.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   Salidas                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V_thres=narrow_vth;
I_drain=i_ds;
I_leak=i_db;
V_early=delta_vth;
end

