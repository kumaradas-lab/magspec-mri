%% This code is a simulation of the T1 Thermometry code with CA (manganese) to try if the correlation between T and T1 works well

t1 = 172.9;
concentration = 1;

% calculate relaxivity r1
t1 = t1*0.001; % convert t1 from ms to s
r1 = 1/t1/concentration; % make relaxivity concentration independant

disp(['Relaxivity r1 = ' num2str(r1) ' Hz']);

%if r1 > 5.38
%    disp('Relaxivity out of range - Unknown parameters');
    
%elseif r1 < 5.38

syms x;
% Multiple equations determined by measurements of differnt
% concentrations all leading to the determination of Temperature
eqn1 = r1 == -0.0013048*x^2 + 0.027533*x + 5.2391; %equation based on 1mM experiments
eqn2 = r1 == -0.00111915*x^2 + 0.021059*x + 5.17005; %equation based on 2mM experiments
eqn3 = r1 == -0.00158477*x^2 + 0.052*x + 4.631; %equation based on 3mM experiments
eqn4 = r1 == -0.0018003*x^2 + 0.0675375*x + 4.4573; %equation based on 4mM experiments
eqn5 = r1 == -0.00130026*x^2 + 0.031408*x + 5.07944; %equation based on 5mM experiments
eqn6 = r1 == -0.00138682*x^2 + 0.03734667*x + 5.03225; %equation based on 6mM experiments

S1 = solve(eqn1, x);
S2 = solve(eqn2, x);
S3 = solve(eqn3, x);
S4 = solve(eqn4, x);
S5 = solve(eqn5, x);
S6 = solve(eqn6, x);

Temperature1 = double(S1);
Temperature2 = double(S2);
Temperature3 = double(S3);
Temperature4 = double(S4);
Temperature5 = double(S5);
Temperature6 = double(S6); % to show S numerically
%Take absolute of temperature for Temperature > 15 deg. C!!
Temperature1 = abs(Temperature1(2));
Temperature2 = abs(Temperature2(2));
Temperature3 = abs(Temperature3(2));
Temperature4 = abs(Temperature4(2));
Temperature5 = abs(Temperature5(2));
Temperature6 = abs(Temperature6(2));

disp(['The Temperature is T1 = ' ,num2str(Temperature1), ' deg. C']);
disp(['The Temperature is T2 = ' ,num2str(Temperature2), ' deg. C']);
disp(['The Temperature is T3 = ' ,num2str(Temperature3), ' deg. C']);
disp(['The Temperature is T4 = ' ,num2str(Temperature4), ' deg. C']);
disp(['The Temperature is T5 = ' ,num2str(Temperature5), ' deg. C']);
disp(['The Temperature is T6 = ' ,num2str(Temperature6), ' deg. C']);

Temperature = (Temperature1+Temperature2+Temperature3+Temperature4+Temperature5+Temperature6)/6;

disp(['The average Temperature is T = ' ,num2str(Temperature), ' deg. C']);


%end

%x = linspace(0, 43, 100)
%y = -0.0013048*x.^2 + 0.027533*x + 5.2391;
%plot(x, y)