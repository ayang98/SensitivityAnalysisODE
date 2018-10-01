function SensitivityAnalysis
Sensitivity(.1,20,2,10,10)

%[t1,xa] = ode15s(@CxGC_SensitivityAnalysis,(0:.1:20),[0 0 0 0 0]);

%plot(t1,xa(:,4))

end


function P_dot = CxGC_SensitivityAnalysis(t,sol) 

%Function containing ChInGON model equations with parameter list for sensitivity analysis

global Pi 


% Concentration of tetracycline and erythromycin
global Tc
global EM 

%Tc = 2; %also use Tc and EM = 10
%EM = 0;

% Induction of the system 
if (t < 200)
    FI_tgene = 1; 
%elseif (t > 250)
        %FI_tgene =1;
else
    FI_tgene = Pi(2);
end

% FIX Parameters
Udil = 0.02; % Cell Dilution
IRES   = 0.5; % IRES contribution

% Protein dynamics ODEs
% P_dot(i,1) = some function of t and sol;

P_dot = zeros(5,1);    

P_dot(1,1) = (Pi(1)*IRES*FI_tgene) - (Pi(15)+Udil)*sol(1); % Total tTA, sol(1)

P_dot(2,1) = (((Pi(4)/(1+(sol(1)/(Pi(8)*(1+(Tc/Pi(13))^Pi(12))))^Pi(12)))... % vhh, sol(2)
             *(1+(Pi(5)*(sol(1)/(Pi(8)*(1+(Tc/Pi(13))^Pi(12))))^Pi(12))))...
             -(Pi(18)+Udil)*sol(2))...
             -(sol(2)*sol(4)*Pi(10))+(sol(5)*Pi(11)); % Proportion of GFP-vhh      

P_dot(3,1) = (((Pi(4)*IRES/(1+(sol(1)/(Pi(8)*(1+(Tc/Pi(13))^Pi(12))))^Pi(12)))... % Total EKRAB, sol(3)
             *(1+(Pi(5)*(sol(1)/(Pi(8)*(1+(Tc/Pi(13))^Pi(12))))^Pi(12))))...
             - (Pi(17)+Udil)*sol(3));...
             
P_dot(4,1) = (((Pi(3)/((1+(sol(1)/(Pi(8)*(1+(Tc/Pi(13))^Pi(12))))^Pi(12))*(1+(sol(3)/(Pi(9)*(1+(EM/Pi(14))))))))... % GFP, sol(4) 
             *(1+(Pi(6)*(sol(1)/(Pi(8)*(1+(Tc/Pi(13))^Pi(12))))^Pi(12))+((Pi(7)*(sol(3)/(Pi(9)*(1+(EM/Pi(14))))))*(1+(sol(1)/(Pi(8)*(1+(Tc/Pi(13))^Pi(12))))^Pi(12)))))...
             -(Pi(16)+Udil)*sol(4))...
             -(sol(2)*sol(4)*Pi(10))+(sol(5)*Pi(11)); % Proportion of GFP-vhh

P_dot(5,1) = (sol(2)*sol(4)*Pi(10))-(sol(5)*Pi(11))-(Pi(19)*Udil)*sol(5); % GFP-vhh complex, sol(5)
         
end

%Note - this function can be used to generate sensitivity coefficient
%graphs for any provided system of ODEs (not just the one above)!
%a is amount of Tc and b is amount of EM
function Sensitivity(delta_t,percentage,display_number, a, b) %display_mode is a number (1 or 2) for number of graphs displayed in a single figure
tic %tic and toc show runtime of the function (displayed in seconds)

names_array = {}; %this is used to label the bar graphs
t0 = 0;
tf = 300;
dt = 1;
t = t0:dt:tf;
sol0 = [0,0,0,0,0]; 
global Pi 
Pi = [10;2;20;60;0.1;10;0.05;3;3;2.8;0.6;2;1;1;0.03;0.03;0.03;0.3;0.3];
global Tc;
global EM;
Tc = a;
EM = b;


bar_graph_final = zeros(length(Pi),2);%this will store all the bar graph y values (every 2 pairs is one group)
[t, sol] = ode15s(@CxGC_SensitivityAnalysis,t,sol0); %recalculate initial conditions
BasalCond = [0,0,0,0,0];  %Basal Conditions 
    for i = 1:length(sol0)
       BasalCond(i) = sol(198,i);
    end 
%disp(BasalCond);

Pi = [10;2;20;60;0.1;10;0.05;3;3;2.8;0.6;2;1;1;0.03;0.03;0.03;0.3;0.3];
Parameter = {'Bo_gene';'FI_gene';'Bo_HP';'Bo_vEK';...
        'FR_TO';'FI_7TO';'FR_EK';'C_tTA';'C_EK';'Kon_vhh';'Koff_vhh';...
            'n';'Kd_Tc';'Kd_EM';'gtTA';'gGFP';'gEK';'gvhh';'gGFPvhh'};

List_Pi = table(Parameter,Pi); % List of parameters with values
for i=1:length(Pi)
    Pi = [10;2;20;60;0.1;10;0.05;3;3;2.8;0.6;2;1;1;0.03;0.03;0.03;0.3;0.3];
    [t1,xa] = ode15s(@CxGC_SensitivityAnalysis,(0:delta_t:400),BasalCond);
    list = zeros(length(t1),2);
    delta_u = (percentage/100)*Pi(i); %delta_u is a percentage of the base P(i) value 
    values_Pi = [Pi(i) Pi(i)+delta_u]; %base value and base value added with delta 
    for j=1:length(values_Pi)
        Pi(i) = values_Pi(j);
        
        [t, sol] = ode15s(@CxGC_SensitivityAnalysis,t,sol0); %recalculate initial conditions depending on value of P(i)
        BasalCond = [0,0,0,0,0];  %Basal Conditions 
        for k = 1:length(sol0)
            BasalCond(k) = sol(198,k);
        end 
        
        %'y(p+delta) and y(p)'disp(Pi)
        [t1,xa] = ode15s(@CxGC_SensitivityAnalysis,(0:delta_t:400), BasalCond);
        %y_list = [y_list xa(i,1)];
        list(:,j) =  xa(:,4); %xa(:,i) where i can range from 1-5 depending on which dependent variable analyzed
    end
    
    bar_graph = [];
    k = 195; %before induction (induction at 200)
    new_list = list(:,1);
    new_list2 = list(:,2);
    for a=1:length(t1)
        if t1(a) == k
            index = a;
            break
        end
    end
    %disp(new_list(index))
    
    while((new_list(index+1) - new_list(index)) > 0.01) %0.001 is the smallest acceptable tolerance!
        index = index+1;
    end
    
    bar_graph =[bar_graph abs((new_list2(index)-new_list(index)))/delta_u];
    
    string1 =   List_Pi{i,{'Parameter'}};
    names_array{end+1} = char(string1);
    
    k2 = 350; %steady state phase after induction  (induction at 200)
    for a2=1:length(t1)
        if t1(a2) == k2
            index2 = a2;
            break
        end
    end
    
    while((new_list(index2+1) - new_list(index2)) > 0.01)
        index2 = index2+1;
    end
    
    bar_graph =[bar_graph abs((new_list2(index2)-new_list(index2)))/delta_u];
    
    bar_graph_final(i,:) = bar_graph;
   
    %disp(new_list(index2))
    final_y = (list(:,2)-list(:,1))/delta_u; 
    if display_number == 2 %display p and p+delta_p in one figure
        figure();
        plot(t1, list(:,1));
        hold on 
        plot(t1, list(:,2));
        string1 = 'GFP(p+delta) and GFP(p) for p = ';
        string2 = List_Pi{i,{'Parameter'}};
        string3 = ' and delta = ';
        string4 = num2str(percentage);
        string5 = 'percent';
        title(strcat(string1, {' '}, string2, {' '}, string3, {' '}, string4, {' '}, string5),'Interpreter', 'none');
        ylabel('y')
        xlabel('time(t)')
        legend('GFP(p)','GFP(p+delta)')
        
    elseif display_number == 1
        figure();
        plot(t1, final_y);
        string3 = 'Sensitivity with respect to ';
        string4 = List_Pi{i,{'Parameter'}};
        string5 = 'with a';
        string6 = num2str(percentage);
        string7 = 'percent increase';
        title(strcat(string3, {' '}, string4, {' '}, string5, {' '}, string6, {' '}, string7),'Interpreter', 'none');
        ylabel('$\frac{d[GFP]}{dp}$','Interpreter','latex','FontSize',14);
        set(get(gca,'ylabel'),'rotation',0)
        xlabel('time(t)')
    end
    
    
end

figure()
%disp(names_array)
c = categorical(names_array);
bar(c,bar_graph_final);
title('Sensitivities of [GFP] to parameters before and after induction')
xlabel('Parameter')
ylabel('d[GFP]/dp')
legend('Before','After')
%set(gca,'YScale','log')
%disp(length(t1))
%disp(length(list(:,1)))

toc
end








