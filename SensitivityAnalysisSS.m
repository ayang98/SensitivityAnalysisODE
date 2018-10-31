function SensitivityAnalysisSS
Sensitivity(20)
%{
global Pi
Pi = [10,2.7,10,500];
[a,b,c]=HM_Optimal
%}
end

%IMPORTANT- right now you have the base parameter optimal triplet as well
%as the delta paramter optimal triplet for each of the 4 paramters (24
%values total) 


function Pvalues = ChInGON_SteadyState(Ind)
global Pi
global Tc
global Em
Parameters = {'Bo_gene';'FI_gene';'B_GFP';'B_vEK'};
% Fold Inductions
FN_tTA = 100;   
FR_tTA  = 0.2;  
FR_EK  = 0.1;   

% Half-life of Proteins 
hf_tTA  = 2;  
hf_EK   = 2;  
hf_vhh  = 2;
hf_GFP  = 24;

% Degradation rates 
gtTA      = 0.693/hf_tTA;
gEK       = 0.693/hf_EK;



% Cell dilution rate 
Udil = 0.02; 

% Trasncription factors binding constants 
KD_tTA  = 1;  
KD_EK   = 1; 

% Dissociation constants [nM]
KD_Tc = 1;       
KD_Em = 1;       

% Influence of IRES 
phi   = 0.5; % FIX - Based on experimental data

% tTA cooperativity 
n_tTA = 2; % FIX - Wenting's paper
n_EK = 2;

if Ind == 0
    FI_tgene = 1;
    gGFP = 0.693/hf_vhh;
else
    FI_tgene = Pi(2);
    gGFP = 0.693/hf_GFP;
end

% Fractions of active transcription factors (tetas)
O_tTA = 1/(1+(Tc/KD_Tc));
O_EK = 1/(1+(Em/KD_Em));

% Rates of production (alphas)

a_tTA = Pi(1)*phi*FI_tgene;
ptTA = a_tTA / (gtTA+Udil); % tTA, sol(1)

a_EK = (Pi(4)/(1+(O_tTA*ptTA/KD_tTA)^n_tTA))*(1+(FR_tTA*(O_tTA*ptTA/KD_tTA)^n_tTA))*phi; 
pEK = a_EK / (gEK+Udil); % EK, sol(2)

a_GFP = (Pi(3)/((1+(O_tTA*ptTA/KD_tTA)^n_tTA)*(1+(O_EK*pEK/KD_EK)^n_EK)))...
            *(1+(FN_tTA*(O_tTA*ptTA/KD_tTA)^n_tTA)+((FR_EK*(O_EK*pEK/KD_EK)^n_EK)*(1+(O_tTA*ptTA/KD_tTA)^n_tTA)));
Pvalues = a_GFP / (gGFP+Udil); % GFP, sol(3)

end



%Note - this function can be used to generate sensitivity coefficient
%graphs for any provided system of ODEs (not just the one above)!
%a is amount of Tc and b is amount of EM
function Sensitivity(percentage) %display_mode is a number (1 or 2) for number of graphs displayed in a single figure
tic %tic and toc show runtime of the function (displayed in seconds);
global Pi
Pi = [10;2.7;10;500];
names_array = {}; %this is used to label the bar graphs
bar_graph_final = zeros(length(Pi),3);
Parameters = {'Bo_gene';'FI_gene';'B_GFP';'B_vEK'};
List_Pi = table(Parameters,Pi); %
for i=1:length(Pi)
    Pi = [10;2.7;10;500];
    delta_u = (percentage/100)*Pi(i); %delta_u is a percentage of the base P(i) value 
    values_Pi = [Pi(i) Pi(i)+delta_u]; %base value and base value added with delta 
    list1 = zeros(1,3);
    list2 = zeros(1,3);
    new_list1 = [];
    new_list2= [];
    for j=1:2
        Pi(i) = values_Pi(j);
        [Maximum_FI, Tc_Opt, Em_Opt] = HM_Optimal;
        if j==1
            list1(1) = Maximum_FI;
            list1(2) = Tc_Opt;
            list1(3) = Em_Opt;
            new_list1 = [new_list1 list1]; %values of triplet with base p
        elseif j==2
            list2(1) = Maximum_FI;
            list2(2) = Tc_Opt;
            list2(3) = Em_Opt; 
            new_list2 = [new_list2 list2]; %values of triplet with delta_p
        end
    end
    string1 =   List_Pi{i,{'Parameters'}};
    names_array{end+1} = char(string1);
    final_list = (new_list2-new_list1)/delta_u;
    bar_graph_final(i,:) = final_list;
    


    
end
%xlswrite('data.xlsx', final_list);
disp(bar_graph_final)
figure()
%disp(names_array)
c = categorical(names_array);
bar(c,bar_graph_final);
title('Sensitivities of the Optimal Tc, EM, and FI to parameters')
xlabel('Parameter')
ylabel('dOpt/dp')
legend('Tc', 'Em', 'FI')
set(gca,'TickLabelInterpreter','none');
%set(gca,'YScale','log')
%disp(length(t1))
%disp(length(list(:,1)))
toc
end


function [Maximum_FI, Tc_Opt, Em_Opt] = HM_Optimal
%When you change the parameters, a new (TC_opt, EM_opt) will be found which
%corresponds to a new FI_opt - we we want to compare the effect changing
%the parameter has on changing the new optimal concentrations and corresponding fold
%indunction
%{
global Bo_gene
global FI_gene
global B_GFP
global B_vEK
%}
global Tc
global Em
%B_GFP = 10;
%B_vEK = 500;
Tc_values_Heat = logspace(0,3,15);   
Tc_values_Plot = logspace(0,2,8);  
Em_values = logspace(1,3,20);   
Bo_values = [10,0.1,3]; %use first value for BIP
FI_values = [2.7,5,2]; %use first value for BIP
tit = {'BIP','ERdj4','EIF4'};
aOut = [length(Tc_values_Heat),length(Em_values)];
ctrl_curve = zeros (length(FI_values),length(Tc_values_Plot));
ind_curve = zeros (length(FI_values),length(Tc_values_Plot));
k = 1; %IMPORTANT - use first value for BIP
%for k = 1:length(FI_values)
%Bo_gene = Bo_values(k);
%FI_gene = FI_values(k);
Maximum_FI = -Inf;
for j = 1:length(Em_values) %loop over all Em values (j is index)
    Em = Em_values(j);
    for i = 1:length(Tc_values_Heat) %loop over all Tc values (i is index)
        Tc = Tc_values_Heat(i);
        for d = [0 1]
            if d == 0
                ctrl = ChInGON_SteadyState(d);
            else 
                induced = ChInGON_SteadyState(d);
            end
        end
        aOut(j,i) = induced/ctrl;
        if aOut(j,i) > Maximum_FI
            Maximum_FI = aOut(j,i);
            Tc_Opt = Tc;
            Em_Opt = Em;
        end
    end
end      
%disp(aOut)
end



