function oval2_problem_new(;
    #u0 = [0.048134401208;0.437983837181;0.00662138171159;0.00536656769503;0.00792812565945;0.0291526890419;0.289047453604;0.000647176292664;0.000147132092929;3.34496998963e-05;7.60461161721e-06;1.72886806243e-06;0.0663154313749;0.908814584851;0.0543600137979;202;297;1.52751139999],
#u0 = [0.128483;1.256853;0.0030203;0.0027977;0.0101511;0.0422942;0.2391346;0.0008014;0.0001464;2.67e-05;4.8e-6;9e-7;0.0619917;1.2444292;0.0486676;199.9383546;137.4267984;1.5180203;1.5180203],
    #u0 = [0.37;0;0;0;0;0;0.21;0;0;0;0;0;0.07;0.03;0.05;4.2;4.2;0],
    u0 = [0.048134401208;0.437983837181;0.00662138171159;0.00536656769503;0.00792812565945;0.0291526890419;0.289047453604;0.000647176292664;0.000147132092929;3.34496998963e-05;7.60461161721e-06;1.72886806243e-06;0.0663154313749;0.908814584851;0.0543600137979;400;297;1.52751139999],

    t_start = 0.0,
    t_final = 5000.0,
    #Parameters

    J1_200=3.,
    J1_34=0.15,
    J2_200=0.2,
    J2_34=0.35,
    J_2z=0.9,#0.75,
    J_O=0.9,
    J0_snail=0.6,
    J1_snail=0.5,
    J2_snail=1.8,
    J1_E=0.1,
    J2_E=0.3,
    J1_V=0.4,
    J2_V=0.4,
    J3_V=2.,
    J1_zeb=3.5,
    J2_zeb=0.9,
    K1=1.0,
    K2=1.0,
    K3=1.0,
    K4=1.0,
    K5=1.0,
    K_TR=20.,
    K_SR=100.,
    TGF0=0.18,#0.2,#0.2,#0.35,#0.1,#0.35,#0.4,#0.3,#0.2,#0.18,#0.18,#0.3,#0.5,#0.4,#0.3,#0.5,#1,#0.18,
    Tk = 1000.,
    k0O=0.35,#0.35,#0.28,#0.35,#0.16,#0.35,
    k0_200=0.0002,
    k0_34=0.001,
    k0_snail=0.0005,
    k0_zeb=0.003,
    k0_TGF=1.1,#1.1,#0.35,#1.1,#1.1,#1.4,#1.3,#1.5,#2,#1.1,
    k0_E=5.,
    k0_V=5.,
    k0_E1=15.,
    k0_E2=5.,
    k0_V1=2.,
    k0_V2=5.,
    k_O=1.2,

    k_200=0.02,
    k_34=0.01,
    k_snail=0.05,
    k_tgf=0.05,
    k_zeb=0.06,
    k_TGF=1.5,
    k_SNAIL=16.,
    k_ZEB=16.,
    kd_ZR1=0.5,
    kd_ZR2=0.5,
    kd_ZR3=0.5,
    kd_ZR4=0.5,
    kd_ZR5=0.5,
    kd_O=1.,
    kd_200=0.035,
    kd_34=0.035,
    kd_SR=0.9,
    kd_E=0.05,
    kd_V=0.05,
    kd_snail=0.09,
    kd_tgf=0.1,
    kd_TR=1.,
    kd_zeb=0.1,
    kd_SNAIL=1.6,
    kd_TGF=0.9,
    kd_ZEB=1.66,


    lamda1=0.5,
    lamda2=0.5,
    lamda3=0.5,
    lamda4=0.5,
    lamda5=0.5,
    lamda_SR=0.5,
    lamda_TR=0.8,

    noise=false,epsilon_m = 0.02ones(18),epsilon_a = zeros(18),
    cell_splitting=false,split_var = 200, split_mean = 500, split_multiple = 0.7)

    d = Dict(:J1_200=>J1_200,
      :u0=>u0,
      :t_start=>t_start,
      :t_final=>t_final,
      :J1_34=>J1_34,
      :J2_200=>J2_200,
      :J2_34=>J2_34,
      :J_2z=>J_2z,
      :J_O=>J_O,
      :J0_snail=>J0_snail,
      :J1_snail=>J1_snail,
      :J2_snail=>J2_snail,
      :J1_E=>J1_E,
      :J2_E=>J2_E,
      #:J_ncad3 => J_ncad3,
      :J1_V=>J1_V,
      :J2_V=>J2_V,
      :J3_V=>J3_V,
      :J1_zeb=>J1_zeb,
      :J2_zeb=>J2_zeb,
      :K1=>K1,
      :K2=>K2,
      :K3=>K3,
      :K4=>K4,
      :K5=>K5,
      :K_TR=>K_TR,
      :K_SR=>K_SR,
      :TGF0=>TGF0,
      #:TGF_flg=>TGF_flg,
      #:Timescale=>Timescale,
      :Tk=>Tk,
      :k0O=>k0O,
      :k0_200=>k0_200,
      :k0_34=>k0_34,
      :k0_snail=>k0_snail,
      :k0_zeb=>k0_zeb,
      :k0_TGF=>k0_TGF,
      :k0_E=>k0_E,
      :k0_V=>k0_V,
      :k0_E1=>k0_E1,
      :k0_E2=>k0_E2,
      :k0_V1=>k0_V1,
      :k0_V2=>k0_V2,
      :k_O=>k_O,
      :k_200=>k_200,
      #:k_SNAIL=>k_SNAIL,
      #:k_TGF=>k_TGF,
      #:k_ZEB=>k_ZEB,
      :k_34=>k_34,
      :k_snail=>k_snail,
      :k_tgf=>k_tgf,
      :k_zeb=>k_zeb,
      :k_TGF=>k_TGF,
      :k_SNAIL=>k_SNAIL,
      :k_ZEB=>k_ZEB,
      :kd_ZR1=>kd_ZR1,
      :kd_ZR2=>kd_ZR2,
      :kd_ZR3=>kd_ZR3,
      :kd_ZR4=>kd_ZR4,
      :kd_ZR5=>kd_ZR5,
      #:kd_SNAIL=>kd_SNAIL,
      :kd_O=>kd_O,
      #:kd_TGF=>kd_TGF,
      #:kd_ZEB=>kd_ZEB,
      :kd_200=>kd_200,
      :kd_34=>kd_34,
      :kd_SR=>kd_SR,
      :kd_E=>kd_E,
      :kd_V=>kd_V,
      :kd_snail=>kd_snail,
      #:kd_Op => kd_Op,
      :kd_tgf=>kd_tgf,
      :kd_zeb=>kd_zeb,
      :kd_SNAIL=>kd_SNAIL,
      :kd_ZEB=>kd_ZEB,
      :kd_TGF=>kd_TGF,
      #:kp_ZEB=>kp_ZEB,
      :lamda1=>lamda1,
      :lamda2=>lamda2,
      :lamda3=>lamda3,
      :lamda4=>lamda4,
      :lamda5=>lamda5,
      :lamda_SR=>lamda_SR,
      :lamda_TR=>lamda_TR)
      #:nO=>nO,
      #:nSO=>nSO,
      #:nzo=>nzo]


    f = function (dy,y,p,t)
      # y[1] = snailt
      # y[2] = SNAIL
      # y[3] = miR34t
      # y[4] = SR1 # abundance of SNAIL/miR34 complex
      # y[5] = zebt
      # y[6] = ZEB
      # y[7] = miR200t
      # y[8] = ZR1 # abundance of ZEB/miR200 complex with i copies of miR200 bound on the sequence of ZEB1
      # y[9] = ZR2
      # y[10] = ZR3
      # y[11] = ZR4
      # y[12] = ZR5
      # y[13] = tgft
      # y[14] = TGF
      # y[15] = tgfR # abundance of TGF/miR200 complex
      # y[16] = Ecad
      # y[17] = Ncad
      # y[18] = Ovol2
      #TGF0=TGF0_coeff*[t>TGF0_cutoff]
      #ODEs
      dy[1]=k0_snail+k_snail*(((y[14]+TGF0)/J0_snail))^2/(1+(((y[14]+TGF0)/J0_snail))^2+(y[18]/J1_snail)^2)/(1+y[2]/J2_snail)-kd_snail*(y[1]-y[4])-kd_SR*y[4];
      dy[2]=k_SNAIL*(y[1]-y[4])-kd_SNAIL*y[2];
      dy[3]=k0_34+k_34/(1+((y[2]/J1_34))^2+((y[6]/J2_34))^2)-kd_34*(y[3]-y[4])-kd_SR*y[4]+lamda_SR*kd_SR*y[4];
      dy[4]=Tk*(K_SR*(y[1]-y[4])*(y[3]-y[4])-y[4]);
      dy[5]=k0_zeb+k_zeb*((y[2]/J1_zeb))^2/(1+((y[2]/J1_zeb))^2+((y[18]/J2_zeb))^6)-kd_zeb*(y[5]-(5*y[8]+10*y[9]+10*y[10]+5*y[11]+y[12]))-kd_ZR1*5*y[8]-kd_ZR2*10*y[9]-kd_ZR3*10*y[10]-kd_ZR4*5*y[11]-kd_ZR5*y[12];
      dy[6]=k_ZEB*(y[5]-(5*y[8]+10*y[9]+10*y[10]+5*y[11]+y[12]))-kd_ZEB*y[6];
      dy[7]=k0_200+k_200/(1+((y[2]/J1_200))^3+((y[6]/J2_200))^2)-kd_200*(y[7]-(5*y[8]+2*10*y[9]+3*10*y[10]+4*5*y[11]+5*y[12])-y[15])-kd_ZR1*5*y[8]-kd_ZR2*2*10*y[9]-kd_ZR3*3*10*y[10]-kd_ZR4*4*5*y[11]-kd_ZR5*5*y[12]+lamda1*kd_ZR1*5*y[8]+lamda2*kd_ZR2*2*10*y[9]+lamda3*kd_ZR3*3*10*y[10]+lamda4*kd_ZR4*4*5*y[11]+lamda5*kd_ZR5*5*y[12]-kd_TR*y[15]+lamda_TR*kd_TR*y[15];
      dy[8]=Tk*(K1*(y[7]-(5*y[8]+2*10*y[9]+3*10*y[10]+4*5*y[11]+5*y[12])-y[15])*(y[5]-(5*y[8]+10*y[9]+10*y[10]+5*y[11]+y[12]))-y[8]);
      dy[9]=Tk*(K2*(y[7]-(5*y[8]+2*10*y[9]+3*10*y[10]+4*5*y[11]+5*y[12])-y[15])*y[8]-y[9]);
      dy[10]=Tk*(K3*(y[7]-(5*y[8]+2*10*y[9]+3*10*y[10]+4*5*y[11]+5*y[12])-y[15])*y[9]-y[10]);
      dy[11]=Tk*(K4*(y[7]-(5*y[8]+2*10*y[9]+3*10*y[10]+4*5*y[11]+5*y[12])-y[15])*y[10]-y[11]);
      dy[12]=Tk*(K5*(y[7]-(5*y[8]+2*10*y[9]+3*10*y[10]+4*5*y[11]+5*y[12])-y[15])*y[11]-y[12]);
      dy[13]=k_tgf-kd_tgf*(y[13]-y[15])-kd_TR*y[15];
      dy[14]=k0_TGF+k_TGF*(y[13]-y[15])-kd_TGF*y[14];
      dy[15]=Tk*(K_TR*(y[7]-(5*y[8]+2*10*y[9]+3*10*y[10]+4*5*y[11]+5*y[12])-y[15])*(y[13]-y[15])-y[15]);
      dy[16]=k0_E+k0_E1/(((y[2]/J1_E))^2+1)+k0_E2/(((y[6]/J2_E))^2+1)-kd_E*y[16];
      dy[17]=k0_V+k0_V1*(((y[2]/J1_V))^2)/(((y[2]/J1_V))^2+1)+k0_V2*(((y[6]/J2_V))^2)/(((y[6]/J2_V))^2+1)/(1+y[18]/J3_V)-kd_V*y[17];
      dy[18]=k0O+k_O/(1+((y[6]/J_O))^2)-kd_O*y[18];

    end
    # sf = StoreFunction(f,d)

    if cell_splitting
      time_choice = (integrator) -> abs(split_var*randn() + split_mean) + integrator.t
      affect! = function (integrator)
        # Apply noise to all
        @. integrator.u .*= max(0.1,1+split_multiple*randn(eltype(integrator.u)))
      end
      cb = IterativeCallback(time_choice,affect!)
    else
      cb = nothing
    end

    if !noise
      return ODEProblem(f,u0,(t_start,t_final),d,callback=cb)
    else
      σ = function (dy,y,p,t)
        dy .= epsilon_m.*y .+ epsilon_a
      end
      return SDEProblem(f,σ,u0,(t_start,t_final),d,callback=cb)
    end
    end

#    struct StoreFunction{F,D} <: Function
#    f::F
#    d::D
#    end
#    (S::StoreFunction)(t,u,du) = S.f(t,u,du)


 #=
 σ = function (t,y,dσ)
  dσ[1] = noiseLevel*1.5y[1]
  dσ[18]= noiseLevel*6y[18]
 end
 =#
