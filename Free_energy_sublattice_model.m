%% Code revision
% This code is developed by Vahid Attari 
% PhD in computational Materials Sciecne
% Affilation: Texas A&M University
% Contact: attari.v@tamu.udu

%% Plot free energy of a phase that is defined by CALPHAD sublattice model


format short

clc
clear all
close all;

% Side fraction ratio
m = 0.25; n = 0.75; R = 8.314;

fig=figure (1);


%% 1st cost function


i = 1;

for T = 1395+273 %:300:1800 %1671
    x_star(:,1:4) = 0;
    for X1 = 0:0.01:1
               
        % Initial guess:
        X0 = [0.18,0.18,0.25,0.45];
    
        ylb = [0,0,0,0];             % lower bound for [y1A,y1B,y2A,y2B]
        yub = [1,1,1,1];             % upper bound for [y1A,y1B,X1,X2]
        
        % linear constraint equations:
        Aeq(1,:) =  [1,0,3,0]; %    [side fraction 1, size fraction 2, side fraction 3, size fraction 4]
        Aeq(2,:) =  [0,1,0,3]; %    [side fraction 1, size fraction 2, side fraction 3, size fraction 4] 
        %Aeq(3,:) = [0.25,0.25,0.75,0.75];   %constraint:m(y1A+y1B)+n(y2C+y2D)=1
        beq = [X1;1-X1];       %[Nb,Ni]
        
        % Cost func: -----> Switch the cost function to see the diffrences
        GTOT = @(x) GTOTfunc(x,R,T);
        %GTOT = @(x) GTOTfunc_2(x,R,T);

        % Composition
        x_star(i,5) = X1;    % Nb
        x_star(i,6) = 1-X1;  % Ni
        options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',5000,'MaxFunctionEvaluations',5000);
        options = optimset('Display', 'off','TolX',1e-8,'TolCon',1e-8) ;
        %x_star(i,1:5) = fmincon(GTOT,X0,[],[],Aeq,beq,ylb,yub,@confun,options);
        [x_star(i,1:4),fval(i),exitflag(i),output(i)] = fmincon(GTOT,X0,[],[],Aeq,beq,ylb,yub,[],options);
        
        i = i+1;
        
    end
    
    for j=1:size(x_star,1)
        GibbsE(j) = GTOTfunc_2(x_star(j,1:4),R,T);
    end
    
    subplot(1,2,1);
    plot(x_star(:,6),GibbsE,'-o','linewidth',1); 
    xlabel('X (Ni)','fontsize',14); ylabel('Gibbs free energy (J/mol)','fontsize',14)
    hold on
    drawnow
end

legend(['NiNb_{3} free energy curve at ',num2str(T),'K'],'fontsize',14);
set(fig,'position',[20 20 800 600])

saveas(fig,'NiNb3_free_energy.jpg')

%% 2nd cost function

% Side fraction ratio
m = 0.25; n = 0.75;  
R = 8.314;

fig=figure (1);

i = 1;

for T = 1395+273 %:300:1800 %1671
    x_star(:,1:4) = 0;
    for X1 = 0:0.01:1
               
        % Initial guess:
        X0 = [0.18,0.18,0.25,0.45];
    
        ylb = [0,0,0,0];             % lower bound for [y1A,y1B,y2A,y2B]
        yub = [1,1,1,1];             % upper bound for [y1A,y1B,X1,X2]
        
        % linear constraint equations:
        Aeq(1,:) =  [1,0,3,0]; %    [side fraction 1, size fraction 2, side fraction 3, size fraction 4]
        Aeq(2,:) =  [0,1,0,3]; %    [side fraction 1, size fraction 2, side fraction 3, size fraction 4] 
        %Aeq(3,:) = [0.25,0.25,0.75,0.75];   %constraint:m(y1A+y1B)+n(y2C+y2D)=1
        beq = [X1;1-X1];       %[Nb,Ni]
        
        % Cost func: -----> Switch the cost function to see the diffrences
        %GTOT = @(x) GTOTfunc(x,R,T);
        GTOT = @(x) GTOTfunc_2(x,R,T);

        % Composition
        x_star(i,5) = X1;    % Nb
        x_star(i,6) = 1-X1;  % Ni
        options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',5000,'MaxFunctionEvaluations',5000);
        options = optimset('Display', 'off','TolX',1e-8,'TolCon',1e-8) ;
        %x_star(i,1:5) = fmincon(GTOT,X0,[],[],Aeq,beq,ylb,yub,@confun,options);
        [x_star(i,1:4),fval(i),exitflag(i),output(i)] = fmincon(GTOT,X0,[],[],Aeq,beq,ylb,yub,[],options);
        
        i = i+1;
        
    end
    
    for j=1:size(x_star,1)
        GibbsE(j) = GTOTfunc_2(x_star(j,1:4),R,T);
    end
    
    subplot(1,2,2);
    plot(x_star(:,6),GibbsE,'-o','linewidth',1); 
    xlabel('X (Ni)','fontsize',14); ylabel('Gibbs free energy (J/mol)','fontsize',14)
    hold on
    drawnow
end

%%

legend(['NiNb_{3} free energy curve at ',num2str(T),'K'],'fontsize',14);
set(fig,'position',[20 20 1100 600])

saveas(fig,'NiNb3_free_energy.jpg')
