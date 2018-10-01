function regression
synthesis(2,20,70)
end

function synthesis(bo,first,second) %bo is the value of the bo_target, first is the first value
%the changing parameter takes, and second is the second value the parameter
%takes
%To-do:
%input bo_target, and the other values are calculated based on experimental
%relationship
%fold induction always 2,5,1.8
%add labels
t0 = 0;
tf = 300;
dt = 1;
t = t0:dt:tf;
sol0 = [0,0,0,0,0]; 
global FN_tTA
FN_tTA = 100;
global FR_tTA
FR_tTA = 0.1;
global FR_EK
FR_EK = 0.5;
global B_GFP
B_GFP = 10;
global B_vEK
%B_vEK = 80;
global Tc
global Em
global Bo_gene
global FI_gene

Tc = 50;
values= [first,second];
for b = 1:length(values)
    B_vEK = values(b);
    array_changes = [bo 2; .15*bo 5; .30*bo 1.8];
    for a = 1:3
        list1 = [];
        list2=[]; 
        Bo_gene = array_changes(a,1);
        FI_gene = array_changes(a,2);
        changing_values = [1 2.5 5 10 20 50 100];
        for i = changing_values
            Em = i; %Tc concentration changing
            [t, sol] = ode15s(@ChInGON_model,t,sol0); %recalculate initial conditions
            BasalCond = [0,0,0,0,0];  %Basal Conditions 
            for k = 1:length(sol0)
               BasalCond(k) = sol(98,k);
            end 
            [t1,xa] = ode15s(@ChInGON_model,(0:.1:1000),BasalCond);
            k1 = 95; %steady state phase before induction  (induction at 100)
            for a1=1:length(t1)
                if t1(a1) == k1
                    index1 = a1;
                    break
                end
            end
            new_list = xa(:,4);
            list1=[list1 new_list(index1)];
            k2 = 800; %steady state phase after induction  (induction at 100)
            for a2=1:length(t1)
                if t1(a2) == k2
                    index2 = a2;
                    break
                end
            end
            list2=[list2 new_list(index2)];
            
        end
        %legend('0','2.5','5','10')
        if b ==1
            k = 0;
        else
            k = 3;
        end
        subplot(2,3,a+k);
        plot(changing_values,list1);
        hold on
        plot(changing_values,list2);
        if ismember(a+k,[1,4])==1
            title(strcat('BIP,',{' '},'Bo_gene=',{' '},num2str(bo)),'Interpreter', 'none')
        elseif ismember(a+k,[2,5])==1
            title('ERdj4')
        elseif ismember(a+k,[3,6])==1
            title('EIF4')
        end
            
        
        
    end
end
end




