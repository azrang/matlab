
% HW 2

close all

%% MAP vs. ML Decision Theory

b1_1 = 1; r1_1 = 9; b2_1 = 1; r2_1 = 9;
case1_MAP = decision_vectors(1, b1_1, r1_1, b2_1, r2_1)
case1_perrMAP = prob_of_error(case1_MAP, b1_1, r1_1, b2_1, r2_1)
case1_ML = decision_vectors(2, b1_1, r1_1, b2_1, r2_1)
case1_perrML = prob_of_error(case1_ML, b1_1, r1_1, b2_1, r2_1)

b1_2 = 4; r1_2 = 6; b2_2 = 1; r2_2 = 9;
case2_MAP = decision_vectors(1, b1_2, r1_2, b2_2, r2_2)
case2_perrMAP = prob_of_error(case2_MAP, b1_2, r1_2, b2_2, r2_2)
case2_ML = decision_vectors(2, b1_2, r1_2, b2_2, r2_2)
case2_perrML = prob_of_error(case2_ML, b1_2, r1_2, b2_2, r2_2)

b1_3 = 1; r1_3 = 9; b2_3 = 4; r2_3 = 6;
case3_MAP = decision_vectors(1, b1_3, r1_3, b2_3, r2_3)
case3_perrMAP = prob_of_error(case3_MAP, b1_3, r1_3, b2_3, r2_3)
case3_ML = decision_vectors(2, b1_3, r1_3, b2_3, r2_3)
case3_perrML = prob_of_error(case3_ML, b1_3, r1_3, b2_3, r2_3)

%% e & f
% The ML decision should not depend on the distribution of the balls in urn 1. 
% This is because the ML decision uses just the likelihood function, which only takes 
% into account the balls in urn 2. This is consistent with the decision rules my 
% algorithm determined, as no matter the distribution of balls between each 
% scenario, the ML decision vector stayed the same.

% The MAP decision should depend on the distribution of balls in Urn 1. Based on 
% the posterior, it takes into account the prior model, which represents the probability
% of red or blue balls in urn 1. This is consistent with what I got from my algorithm, 
% as MAP differed from case to case. 

%% simulation
%b1 is 1, r1 is 2, b2 is 3, r2 is 4
x1_b1= [3 3 4 4 4 4 4 4 4 4 4]; %choose b1 and place in urn2
x1_r1= [3 4 4 4 4 4 4 4 4 4 4]; %choose r1 and place in urn2
x2_b1 = [3 3 4 4 4 4 4 4 4 4 4 ]; %choose b1 and place in urn2
x2_r1 = [3 4 4 4 4 4 4 4 4 4 4 ]; 
x3_b1 = [3 3 3 3 3 4 4 4 4 4 4];
x3_r1 = [3 3 3 3 4 4 4 4 4 4 4];   %URN2

err_MAP_case1=0; %MAP
err_MAP_case2=0;
err_MAP_case3=0;
err_ML_case1=0; %ML index
err_ML_case2=0;
err_ML_case3=0;

% randomly choose which color from urn 1, if blue from urn 1, then 
% choose from x1_b1

exp_decision_1 =[];
exp_decision_2 =[];
exp_decision_3 =[];

for N = 1:10e5
        ch1_b1 = randi(length(x1_b1));  %choosing from urn 2
        choi1_b1 = x1_b1(ch1_b1);
        ch1_r1 = randi(length(x1_r1));
        choi1_r1 = x1_r1(ch1_r1);   %if value is 2
        if(choi1_b1 == 3) 
            exp_decision_1(1) = 1; %blue from urn1
        elseif(choi1_b1 == 4) %red
            exp_decision_1(2) =1; %blue from urn1
        end
        if(choi1_r1 == 3)
            exp_decision_1(1) = 2;
        elseif(choi1_b1 == 4) %red
            exp_decision_1(2) =2;
        end

        ch2_b1 = randi(length(x2_b1));
        choi2_b1 = x2_b1(ch2_b1);
        ch2_r1 = randi(length(x2_r1));
        choi2_r1 = x2_r1(ch2_r1);   %if value is 2
        if(choi2_b1 == 3) 
            exp_decision_2(1) = 1;
        elseif(choi2_b1 == 4) %r2
            exp_decision_2(2) = 1;
        end
        if(choi2_r1 == 3) %b2
            exp_decision_2(1) = 2;  %choose red
        elseif(choi2_r1 == 4) %r2
            exp_decision_2(2) = 2;
        end

        ch3_b1 = randi(length(x3_b1));
        choi3_b1 = x3_b1(ch3_b1);
        ch3_r1 = randi(length(x3_r1));
        choi3_r1 = x3_r1(ch3_r1);   %if value is 2
        if(choi3_b1 == 3) %b2
            exp_decision_3(1) = 1;
        elseif(choi3_b1 == 4) %r2
            exp_decision_3(2) = 1;
        end
        if(choi3_r1 == 3) %b2
            exp_decision_3(1) = 2;  %choose red
        elseif(choi3_r1 == 4) %r2
            exp_decision_3(2) = 2;
        end

     if(exp_decision_1(1) ~= 2 && exp_decision_1(2) ~= 2)   %MAP [22] ML [12]
        err_MAP_case1 = err_MAP_case1 +1;
     end
     if(exp_decision_1(1) ~= 1 && exp_decision_1(2) ~= 2 )
         err_ML_case1 = err_ML_case1 +1;
     end

     if(exp_decision_2(1) ~= 1 && exp_decision_2(2) ~= 2)%[1 2]
         err_MAP_case2 = err_MAP_case2 +1;
         err_ML_case2 =err_ML_case2 +1;
     end

   if(exp_decision_3(1) ~= 1 && exp_decision_3(2) ~= 2)%MAP[2 2], ML[12]
         err_ML_case3 =err_ML_case3 +1;
   elseif(exp_decision_3(1) ~= 2 && exp_decision_3(2) ~= 2)
        err_MAP_case3 = err_MAP_case3 +1;
   end
end

disp('ML errors case 1')
disp(err_ML_case1/10e5)
disp('ML errors case 2')
disp(err_ML_case2/10e5)
disp('ML errors case 3')
disp(err_ML_case3/10e5)
disp('MAP errors case 1')
disp(err_MAP_case1/10e5)
disp('MAP errors case 2')
disp(err_MAP_case2/10e5)
disp('MAP errors case 3')
disp(err_MAP_case3/10e5)

% the experimental errors found for both MAP and ML were lower than the 
% theoretically computed errors. For example, for case 1, ML experiemntal 
%error was .1 lower than theoretical and the MAP was also around .1 lower. 
% This pattern was similar for all the cases. For case 2, the experimental
% errors were much lower than the theoretical. 

%% functions 

function [pi_r, pi_b] = compute_prior(b1, r1)
    pi_r = r1/(r1+b1);
    pi_b = b1/(r1+b1);
end 

function [lr2r1, lb2b1, lr2b1, lb2r1]= compute_likelihood(b2, r2)  %number of balls in each urn
    lr2r1 = (r2+1)/ (r2+b2+1);
    lb2r1 = b2/(r2+b2+1); 
    lr2b1 = r2/(b2+r2+1);
    lb2b1 = (b2+1)/(b2+1+r2);
end

function [posr1r2, posb1b2, posr1b2, posb1r2]= compute_posteriori(b1, r1, b2, r2) 
    pi_r = r1/(r1+b1);
    pi_b = b1/(r1+b1);
    posr1r2 = (((r2+1)/ (r2+b2+1)) * pi_r) / (((r2+1)/(r2+b2+1)) * pi_r + (r2/(b2+r2+1))* pi_b );
    posb1b2 = (((b2+1)/(b2+1+r2)) * pi_b) / (((b2+1)/(b2+1+r2)) * pi_b + (b2/(r2+b2+1)) *pi_r) ;
    posr1b2 = ((b2/(r2+b2+1)) * pi_r) /  ((b2/(r2+b2+1)) * pi_r + ((b2+1)/(b2+1+r2)) *pi_b) ;
    posb1r2 = ((r2/(b2+r2+1)) * pi_b) / ((r2/(b2+r2+1)) * pi_b + ((r2+1)/ (r2+b2+1)) * pi_r) ;
end
    
function v= decision_vectors(type, b1, r1, b2, r2)  %type determines MAP or ML
    lr2r1 = (r2+1)/ (r2+b2+1);
    lb2r1 = b2/(r2+b2+1); 
    lr2b1 = r2/(b2+r2+1);
    lb2b1 = (b2+1)/(b2+1+r2);
    pi_r = r1/(r1+b1);
    pi_b = b1/(r1+b1);
    posr1r2 = (((r2+1)/ (r2+b2+1)) * pi_r) / (((r2+1)/(r2+b2+1)) * pi_r + (r2/(b2+r2+1))* pi_b );
    posb1b2 = (((b2+1)/(b2+1+r2)) * pi_b) / (((b2+1)/(b2+1+r2)) * pi_b + (b2/(r2+b2+1)) *pi_r) ;
    posr1b2 = ((b2/(r2+b2+1)) * pi_r) /  ((b2/(r2+b2+1)) * pi_r + ((b2+1)/(b2+1+r2)) *pi_b) ;
    posb1r2 = ((r2/(b2+r2+1)) * pi_b) / ((r2/(b2+r2+1)) * pi_b + ((r2+1)/ (r2+b2+1)) * pi_r) ;

    if type == 1  %MAP     P(B1|B2)= P(B2|B1)P(B1)  vs P(R1|B2)= P(B2|R1)P(R1)
        if posb1b2 > posr1b2
            v(1) = 1;  %choosing what decision from urn 1 (BLUE)
        else 
            v(1) = 2;  %choosing RED
        end
        if posb1r2 > posr1r2   %P(R1|R2) vs P(B1|R2)
            v(2) = 1;
        else 
            v(2) = 2;
        end
    elseif type==2  %ML, max OF LIKELIHOOD
        if (lb2b1) > (lb2r1)
            v(1) = 1;
        elseif (lb2b1) < (lb2r1)    %P(urn2 | urn1)
            v(1) = 2;
        end
        if (lr2b1) > (lr2r1)
            v(2) = 1;
        elseif (lr2b1) < (lr2r1)
            v(2) = 2; 
        end
    end
end

function pe = prob_of_error(v, b1, r1, b2, r2)
    lr2r1 = (r2+1)/ (r2+b2+1);
    lb2r1 = b2/(r2+b2+1); 
    lr2b1 = r2/(b2+r2+1);
    lb2b1 = (b2+1)/(b2+1+r2);
    pi_r = r1/(r1+b1);
    pi_b = b1/(r1+b1);
    if v(1) == 1  %B2 | B1
        if v(2) == 1 %R2 | B1
            %B2&R1 + R2&R1
            pe = lb2r1 * pi_r + lr2r1 * pi_r;
        elseif v(2) == 2 %R2 | R1
            %B2&R1 + R2&B1
            pe = lb2r1*pi_r + lr2b1*pi_b;
        end
    elseif v(1) == 2 %B2 | R1
        if v(2) == 1 %R2 | B1
            %B2&B1 + R2&R1
            pe = lb2b1 * pi_b + lr2r1*pi_r;
        elseif v(2) == 2 %R2 | R1
            %B2&B1 + R2&B1
            pe = lb2b1 *pi_b + lr2b1*pi_b;
        end
    end
end


