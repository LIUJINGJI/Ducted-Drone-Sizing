%%%%%%%%%momentum approach for duct rotor design%%%%%%%%%%%%%
%%%%%%%%%%%%%Jingjiang Liu @ HKU 25-09-2024%%%%%%%%%%%%%%%%%%
%V_0 = 0; %Forward Velocity
%V_1 = ; %Disk output Velocity
clc
clear
%%Parameter define 
k_rotor = 1; %thurst ratio
k_fan = 1 - k_rotor; %thurst ratio
K_1 = 0; %loss of intake
K_2 = 1; %loss of kinetic

error = 0.005;% define the minimal error

blades_number = 3; %blades_number动量法中没有意义

rpm = 6000; %rotation speed

height = 0.10; %designed height of duct

skin = 1.2 / 1000;

materia_denisty = 1200;

thick_factor = 1.2;%初始1.2

Weight_ratio0 = 10;%迭代精度初始化

k = 1;

enduranceTime=0.1; %设计续航时间h

eff=6;%力效g/w,越高越好

weightpayload = 0.3; %载荷重量kg

otherweightfactor = 0.1; %其他项目占比


for weight_guess = 0.3:0.01:2

Thrust_P = 2*weight_guess * 9.8; %Designed Thrust (N)


weightbattery = ((weight_guess*1000/eff) * enduranceTime)/150;

%Duct and Rotor
%D_rotor = 0.2; %Disk dia (m)
i= 1;
for D_rotor = 0.05:0.001:0.4 %search space and step size

R_rotor = D_rotor/2; %Blade length (m)

A_1 = pi * R_rotor^2; %Blade Disk Area (m/s)

A_2 = 1.0 * A_1; %Duct Output Area

rho = 1.25; %sealevel denisity (kg/m³)

%Fluid Velocity
V_0 = 0; %Forward Velocity

V_2 = sqrt(Thrust_P/(rho*A_2)); %Duct output Velocity

Q = A_2 * V_2; %volume of air

V_1 = Q/A_1; %Velocity at Disk

p_1_delta = k_rotor * Thrust_P / A_1; %statistic pressure differential loss excluded

p_delta = p_1_delta + 0.5 * rho * (V_1^2); %statistic pressure differential loss inxcluded

P_tF = p_delta + 0.5 * rho * (V_1^2); %Pressure total

n_s = rpm * ((Q^0.5) / (P_tF^0.75)); %rpm scale

K_u = 0.0060715 * n_s + 0.97143; %A ration parameter no meaning

D = (60 * K_u * sqrt(2*P_tF/rho))/(pi * rpm);


if D_rotor - D  < error

   D_good = D_rotor;

   D_good2 = D;

   V_FIANL = V_1;

end

convergenceRecorder(k,1)= D_good - D_rotor;
   
end

convergence = D_good - D_rotor;

Weight_stru = D_good * pi * height * thick_factor * 2 * skin * materia_denisty;

maxpower = 2*weight_guess*1000/eff;

powerweight = (-0.0002*maxpower^2 + 0.3982*maxpower - 12.876)/1000;

Weight_take_off = (Weight_stru + weightbattery + weightpayload+powerweight)/(1-otherweightfactor);

Weight_ratio = abs(Weight_take_off - weight_guess)/Weight_take_off;

Weight_ratio_recorder(k,1) = Weight_ratio;

Weight_ratio_recorder(k,2) = weight_guess;

Weight_ratio_recorder(k,3) = Weight_take_off;

Weight_ratio_recorder(k,4) = weightbattery;

Weight_ratio_recorder(k,5) = D_good2;

Weight_ratio_recorder(k,6) = Weight_take_off-weight_guess;

if Weight_ratio < Weight_ratio0

    final_weight = Weight_take_off; 

    final_weight_guess = weight_guess;

    Weight_ratio0 =  Weight_ratio; %renew the ratio

    final_D_good = D_good;

    finalBatt = weightbattery;

    finalpowerweight = powerweight;

    finalmaxpower = maxpower;
end

k = k +1;
end

final_weight_guess

final_D_good

finalBatt

finalpowerweight

finalmaxpower

%% 绘制收敛图像（Weight Guess vs Weight Ratio Error）
figure;
plot(Weight_ratio_recorder(:,2), Weight_ratio_recorder(:,1), 'o-', 'LineWidth', 1.5);
grid on;
xlabel('Weight Guess (kg)', 'FontSize', 12);
ylabel('Weight Ratio Error', 'FontSize', 12);
title('Iteration Convergence of Takeoff Weight Estimation', 'FontSize', 14);
%% 绘制重量估算误差收敛图像
figure;
plot(Weight_ratio_recorder(:,2), Weight_ratio_recorder(:,6), '^-', ...
    'LineWidth', 1.5, 'Color', [0.85 0.33 0.1]);
grid on;

xlabel('Weight Guess (kg)', 'FontSize', 12);
ylabel('Weight Error (kg)', 'FontSize', 12);
title('Convergence of Takeoff Weight Estimation', 'FontSize', 14);
set(gca, 'FontSize', 11);
