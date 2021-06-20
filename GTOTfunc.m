function GTOT = GTOTfunc(y,R,T)

  
   %%%%%% Free energy of NbNi3 ----> Source unknown
   
   y1A = y(1);
   y1B = y(2);  
   y2C = y(3);
   y2D = y(4);  
   
   m = 0.25;
   n = 0.75;
   
   H0_Nb_beta = 0;
   H0_Ni_gamma= 0;
   G0_Nb_beta = -8519.353 + 142.045475*T - 26.4711*T*log(T) + 2.03475e-4*T^2 - 3.5012e-7*T^3 + 93399/T + H0_Nb_beta; %valid between 298 and 2750K
   G0_Ni_gamma= -5179.159 + 117.854*T - 22.096*T*log(T) - 0.0048407*T^2 + H0_Ni_gamma;                               %valid between 298 and 1728K
   
   G0_Ni_beta = +8715.074 - 3.556*T + G0_Ni_gamma;
   G0_Nb_gamma= +13500 + 1.7*T + G0_Nb_beta;
   
   % 
   
   G0_NbN3_NbNi = 0.25*G0_Nb_beta + 0.75*G0_Ni_gamma - 35300.600*T + 4.83322*T;
   G0_NbN3_NiNb = 0.25*G0_Ni_gamma+ 0.75*G0_Nb_beta  + 45300.575*T - 4.83322*T;
   G0_NbN3_NbNb = G0_Nb_beta  + 5000;
   G0_NbN3_NiNi = G0_Ni_gamma + 5000;
   
   L_NbNi3_NbNiNb = -3079.625;
   L_NbNi3_NbNiNi = -3079.625;   % L_NbNi3_NbNiNi
   L_NbNi3_NbNbNi = +13505.625;  % L_NbNi3_NbNbNi
   L_NbNi3_NiNbNi = +13505.625;
   
   F_R = y1A.*y2C.*G0_NbN3_NbNb + y1B.*y2C.*G0_NbN3_NiNb + y1A.*y2D.*G0_NbN3_NbNi + y1B.*y2D.*G0_NbN3_NiNi;
   F_M = m*R*T*( y1A.*log(y1A) + y1B.*log(y1B) ) + n*R*T*( y2C.*log(y2C) + y2D.*log(y2D) );
   F_ex= y1A.*y1B.*y2C.*(L_NbNi3_NbNiNb) + y1A.*y1B.*y2D.*(L_NbNi3_NbNiNi) + y2C.*y2D.*y1A.*(L_NbNi3_NbNbNi) + y2C.*y2D.*y1B.*(L_NbNi3_NiNbNi);
   GTOT= F_R + F_M + F_ex;

   
end

