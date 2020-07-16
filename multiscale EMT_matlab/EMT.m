function [yat] = EMT(u0,TGF0,k0O,t_start,t_final,epsilon_m,saveat)
    dt=0.0001; %step size
    pd = makedist('Normal',0,sqrt(dt)); % Initialize the probability distribution for our 
                         % random variable with mean 0 and 
                         % stdev of sqrt(dt)
   
    
    N=(t_final-t_start)/dt+1;

    ts = linspace(t_start, t_final, N);
%     nc = size(u0,1); % number of cells
%     y = zeros(nc,18,N);
%     y(:,:,1)=u0;
%     for i = 2:numel(ts)
%         %t = (i-1) .* dt;
%         dW  = random(pd);
%         y(:,:,i) = y(:,:,i-1) + emtdy(y(:,:,i-1),TGF0,k0O) .* dt + epsilon_m.*y(:,:,i-1).* dW;
%     end
    y=u0;
    j=1;
    for i = 2:numel(ts)
        t = (i-1) .* dt;
        dW  = random(pd);
        y =  max(0, y + emtdy(y,TGF0,k0O) .* dt + epsilon_m.*y.* dW);
        if mod(t,saveat)==0
            yat(:,:,j) = y;
            j=j+1;
        end
            
    end
    


    
    
    function dy = emtdy(y,TGF0,k0O)
        if isempty(TGF0)
            TGF0 = 0.35;
        end
        if isempty(k0O)
            k0O = 0.35;
        end

        J1_200=3.;
        J1_34=0.15;
        J2_200=0.2;
        J2_34=0.35;
        J_2z=0.9;%0.75,
        J_O=0.9;
        J0_snail=0.6;
        J1_snail=0.5;
        J2_snail=1.8;
        J1_E=0.1;
        J2_E=0.3;
        J1_V=0.4;
        J2_V=0.4;
        J3_V=2.;
        J1_zeb=3.5;
        J2_zeb=0.9;
        K1=1.0;
        K2=1.0;
        K3=1.0;
        K4=1.0;
        K5=1.0;
        K_TR=20.;
        K_SR=100.;
        %TGF0=0.35; %#0.2,#0.2,#0.35,#0.1,#0.35,#0.4,#0.3,#0.2,#0.18,#0.18,#0.3,#0.5,#0.4,#0.3,#0.5,#1,#0.18,
        Tk = 1000.;
        %k0O=0.35;%#0.35,#0.28,#0.35,#0.16,#0.35,
        k0_200=0.0002;
        k0_34=0.001;
        k0_snail=0.0005;
        k0_zeb=0.003;
        k0_TGF=1.1;%#1.1,#0.35,#1.1,#1.1,#1.4,#1.3,#1.5,#2,#1.1,
        k0_E=5.;
        k0_V=5.;
        k0_E1=15.;
        k0_E2=5.;
        k0_V1=2.;
        k0_V2=5.;
        k_O=1.2;

        k_200=0.02;
        k_34=0.01;
        k_snail=0.05;
        k_tgf=0.05;
        k_zeb=0.06;
        k_TGF=1.5;
        k_SNAIL=16.;
        k_ZEB=16.;
        kd_ZR1=0.5;
        kd_ZR2=0.5;
        kd_ZR3=0.5;
        kd_ZR4=0.5;
        kd_ZR5=0.5;
        kd_O=1.;
        kd_200=0.035;
        kd_34=0.035;
        kd_SR=0.9;
        kd_E=0.05;
        kd_V=0.05;
        kd_snail=0.09;
        kd_tgf=0.1;
        kd_TR=1.;
        kd_zeb=0.1;
        kd_SNAIL=1.6;
        kd_TGF=0.9;
        kd_ZEB=1.66;


        lamda1=0.5;
        lamda2=0.5;
        lamda3=0.5;
        lamda4=0.5;
        lamda5=0.5;
        lamda_SR=0.5;
        lamda_TR=0.8;


          % y[1] = snailt
          % y[2] = SNAIL
          % y[3] = miR34t
          % y[4] = SR1 # abundance of SNAIL/miR34 complex
          % y[5] = zebt
          % y[6] = ZEB
          % y[7] = miR200t
          % y[8] = ZR1 # abundance of ZEB/miR200 complex with i copies of miR200 bound on the sequence of ZEB1
          % y[9] = ZR2
          % y[10] = ZR3
          % y[11] = ZR4
          % y[12] = ZR5
          % y[13] = tgft
          % y[14] = TGF
          % y[15] = tgfR # abundance of TGF/miR200 complex
          % y[16] = Ecad
          % y[17] = Ncad
          % y[18] = Ovol2


%         dy(1)=k0_snail+k_snail*(((y(14)+TGF0)/J0_snail))^2/(1+(((y(14)+TGF0)/J0_snail))^2+(y(18)/J1_snail)^2)/(1+y(2)/J2_snail)-kd_snail*(y(1)-y(4))-kd_SR*y(4);
%         dy(2)=k_SNAIL*(y(1)-y(4))-kd_SNAIL*y(2);
%         dy(3)=k0_34+k_34/(1+((y(2)/J1_34))^2+((y(6)/J2_34))^2)-kd_34*(y(3)-y(4))-kd_SR*y(4)+lamda_SR*kd_SR*y(4);
%         dy(4)=Tk*(K_SR*(y(1)-y(4))*(y(3)-y(4))-y(4));
%         dy(5)=k0_zeb+k_zeb*((y(2)/J1_zeb))^2/(1+((y(2)/J1_zeb))^2+((y(18)/J2_zeb))^6)-kd_zeb*(y(5)-(5*y(8)+10*y(9)+10*y(10)+5*y(11)+y(12)))-kd_ZR1*5*y(8)-kd_ZR2*10*y(9)-kd_ZR3*10*y(10)-kd_ZR4*5*y(11)-kd_ZR5*y(12);
%         dy(6)=k_ZEB*(y(5)-(5*y(8)+10*y(9)+10*y(10)+5*y(11)+y(12)))-kd_ZEB*y(6);
%         dy(7)=k0_200+k_200/(1+((y(2)/J1_200))^3+((y(6)/J2_200))^2)-kd_200*(y(7)-(5*y(8)+2*10*y(9)+3*10*y(10)+4*5*y(11)+5*y(12))-y(15))-kd_ZR1*5*y(8)-kd_ZR2*2*10*y(9)-kd_ZR3*3*10*y(10)-kd_ZR4*4*5*y(11)-kd_ZR5*5*y(12)+lamda1*kd_ZR1*5*y(8)+lamda2*kd_ZR2*2*10*y(9)+lamda3*kd_ZR3*3*10*y(10)+lamda4*kd_ZR4*4*5*y(11)+lamda5*kd_ZR5*5*y(12)-kd_TR*y(15)+lamda_TR*kd_TR*y(15);
%         dy(8)=Tk*(K1*(y(7)-(5*y(8)+2*10*y(9)+3*10*y(10)+4*5*y(11)+5*y(12))-y(15))*(y(5)-(5*y(8)+10*y(9)+10*y(10)+5*y(11)+y(12)))-y(8));
%         dy(9)=Tk*(K2*(y(7)-(5*y(8)+2*10*y(9)+3*10*y(10)+4*5*y(11)+5*y(12))-y(15))*y(8)-y(9));
%         dy(10)=Tk*(K3*(y(7)-(5*y(8)+2*10*y(9)+3*10*y(10)+4*5*y(11)+5*y(12))-y(15))*y(9)-y(10));
%         dy(11)=Tk*(K4*(y(7)-(5*y(8)+2*10*y(9)+3*10*y(10)+4*5*y(11)+5*y(12))-y(15))*y(10)-y(11));
%         dy(12)=Tk*(K5*(y(7)-(5*y(8)+2*10*y(9)+3*10*y(10)+4*5*y(11)+5*y(12))-y(15))*y(11)-y(12));
%         dy(13)=k_tgf-kd_tgf*(y(13)-y(15))-kd_TR*y(15);
%         dy(14)=k0_TGF+k_TGF*(y(13)-y(15))-kd_TGF*y(14);
%         dy(15)=Tk*(K_TR*(y(7)-(5*y(8)+2*10*y(9)+3*10*y(10)+4*5*y(11)+5*y(12))-y(15))*(y(13)-y(15))-y(15));
%         dy(16)=k0_E+k0_E1/(((y(2)/J1_E))^2+1)+k0_E2/(((y(6)/J2_E))^2+1)-kd_E*y(16);
%         dy(17)=k0_V+k0_V1*(((y(2)/J1_V))^2)/(((y(2)/J1_V))^2+1)+k0_V2*(((y(6)/J2_V))^2)/(((y(6)/J2_V))^2+1)/(1+y(18)/J3_V)-kd_V*y(17);
%         dy(18)=k0O+k_O/(1+((y(6)/J_O))^2)-kd_O*y(18);
        dy(:,1)=k0_snail+k_snail*(((y(:,14)+TGF0)/J0_snail)).^2./(1+(((y(:,14)+TGF0)/J0_snail)).^2+(y(:,18)/J1_snail).^2)./(1+y(:,2)/J2_snail)-kd_snail*(y(:,1)-y(:,4))-kd_SR*y(:,4);
        dy(:,2)=k_SNAIL*(y(:,1)-y(:,4))-kd_SNAIL*y(:,2);
        dy(:,3)=k0_34+k_34./(1+((y(:,2)/J1_34)).^2+((y(:,6)/J2_34)).^2)-kd_34*(y(:,3)-y(:,4))-kd_SR*y(:,4)+lamda_SR*kd_SR*y(:,4);
        dy(:,4)=Tk*(K_SR*(y(:,1)-y(:,4)).*(y(:,3)-y(:,4))-y(:,4));
        dy(:,5)=k0_zeb+k_zeb*((y(:,2)/J1_zeb)).^2./(1+((y(:,2)/J1_zeb)).^2+((y(:,18)/J2_zeb)).^6)-kd_zeb*(y(:,5)-(5*y(:,8)+10*y(:,9)+10*y(:,10)+5*y(:,11)+y(:,12)))-kd_ZR1*5*y(:,8)-kd_ZR2*10*y(:,9)-kd_ZR3*10*y(:,10)-kd_ZR4*5*y(:,11)-kd_ZR5*y(:,12);
        dy(:,6)=k_ZEB*(y(:,5)-(5*y(:,8)+10*y(:,9)+10*y(:,10)+5*y(:,11)+y(:,12)))-kd_ZEB*y(:,6);
        dy(:,7)=k0_200+k_200./(1+((y(:,2)/J1_200)).^3+((y(:,6)/J2_200)).^2)-kd_200*(y(:,7)-(5*y(:,8)+2*10*y(:,9)+3*10*y(:,10)+4*5*y(:,11)+5*y(:,12))-y(:,15))-kd_ZR1*5*y(:,8)-kd_ZR2*2*10*y(:,9)-kd_ZR3*3*10*y(:,10)-kd_ZR4*4*5*y(:,11)-kd_ZR5*5*y(:,12)+lamda1*kd_ZR1*5*y(:,8)+lamda2*kd_ZR2*2*10*y(:,9)+lamda3*kd_ZR3*3*10*y(:,10)+lamda4*kd_ZR4*4*5*y(:,11)+lamda5*kd_ZR5*5*y(:,12)-kd_TR*y(:,15)+lamda_TR*kd_TR*y(:,15);
        dy(:,8)=Tk*(K1*(y(:,7)-(5*y(:,8)+2*10*y(:,9)+3*10*y(:,10)+4*5*y(:,11)+5*y(:,12))-y(:,15)).*(y(:,5)-(5*y(:,8)+10*y(:,9)+10*y(:,10)+5*y(:,11)+y(:,12)))-y(:,8));
        dy(:,9)=Tk*(K2*(y(:,7)-(5*y(:,8)+2*10*y(:,9)+3*10*y(:,10)+4*5*y(:,11)+5*y(:,12))-y(:,15)).*y(:,8)-y(:,9));
        dy(:,10)=Tk*(K3*(y(:,7)-(5*y(:,8)+2*10*y(:,9)+3*10*y(:,10)+4*5*y(:,11)+5*y(:,12))-y(:,15)).*y(:,9)-y(:,10));
        dy(:,11)=Tk*(K4*(y(:,7)-(5*y(:,8)+2*10*y(:,9)+3*10*y(:,10)+4*5*y(:,11)+5*y(:,12))-y(:,15)).*y(:,10)-y(:,11));
        dy(:,12)=Tk*(K5*(y(:,7)-(5*y(:,8)+2*10*y(:,9)+3*10*y(:,10)+4*5*y(:,11)+5*y(:,12))-y(:,15)).*y(:,11)-y(:,12));
        dy(:,13)=k_tgf-kd_tgf*(y(:,13)-y(:,15))-kd_TR*y(:,15);
        dy(:,14)=k0_TGF+k_TGF*(y(:,13)-y(:,15))-kd_TGF*y(:,14);
        dy(:,15)=Tk*(K_TR*(y(:,7)-(5*y(:,8)+2*10*y(:,9)+3*10*y(:,10)+4*5*y(:,11)+5*y(:,12))-y(:,15)).*(y(:,13)-y(:,15))-y(:,15));
        dy(:,16)=k0_E+k0_E1./(((y(:,2)/J1_E)).^2+1)+k0_E2./(((y(:,6)/J2_E)).^2+1)-kd_E*y(:,16);
        dy(:,17)=k0_V+k0_V1*(((y(:,2)/J1_V)).^2)./(((y(:,2)/J1_V)).^2+1)+k0_V2*(((y(:,6)/J2_V)).^2)./(((y(:,6)/J2_V)).^2+1)./(1+y(:,18)/J3_V)-kd_V*y(:,17);
        dy(:,18)=k0O+k_O./(1+((y(:,6)/J_O)).^2)-kd_O*y(:,18);
    end

end

