function GTOT = GTOTfunc_2(y,R,T)

   %%%%%% Data taken from Assessment of the Nb-Ni system by joubert et al (2004)
   
   yNb_u = y(1);
   yNi_u = y(2);  
   yNb_v = y(3);
   yNi_v = y(4);
   
   G_Nb_bcc = -8519.353 + 142.045475*T - 26.4711*T*log(T) + 0.203475e-3*T^2 - 0.35012e-6*T^3 + 93399*T^(-1);         
   G_Ni_fcc = -5179.159 + 117.854*T - 22.096*T*log(T) - 4.8407E-3*T^2; %+ Gpres + Gmag; %(298.15 < T < 1728)                                                                                  
   
   a_NbNi3_NbNi = -128000; b_NbNi3_NbNi = 33;

   G_NbNi3_NbNi = a_NbNi3_NbNi + b_NbNi3_NbNi*T + G_Nb_bcc + 3*G_Ni_fcc;
   G_NbNi3_NbNb = 20000 + 4*G_Nb_bcc;
   G_NbNi3_NiNi = 20000 + 4*G_Ni_fcc;
   G_NbNi3_NiNb = 40000 - a_NbNi3_NbNi - b_NbNi3_NbNi*T + G_Ni_fcc + 3*G_Nb_bcc;
   
   L_NbNi3_NbNi_S = 100000;
   L_NbNi3_S_NbNi = 0;
   
   G_ref_NbNi3  = yNb_u*yNb_v*G_NbNi3_NbNb + yNi_u*yNi_v*G_NbNi3_NiNi + yNb_u*yNi_v*G_NbNi3_NbNi + yNi_u*yNb_v*G_NbNi3_NiNb ;
   G_ide_NbNi3  = R*T*( yNb_u.*log(yNb_u) + yNi_u.*log(yNi_u) ) + 3*R*T*( yNb_v.*log(yNb_v) + yNi_v.*log(yNi_v) );
   G_exe_NbNi3  = yNb_u*yNi_u*L_NbNi3_NbNi_S + yNb_v*yNi_v*L_NbNi3_S_NbNi;
   
   GTOT = G_ref_NbNi3 + G_ide_NbNi3 + G_exe_NbNi3;
   
   
end 