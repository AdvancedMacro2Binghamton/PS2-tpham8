%%%Setting up parameters%%%
delta= 0.025;
alpha= 0.35;
beta= 0.99;
sigma= 2;
tmtx= [0.977 0.023; 0.074 0.926];
A=[1.1 0.678];
num_k=100
%%%State Space%%%
k_min= 0;
k_max= 45
k= linspace(k_min, k_max, num_k); 

kmat = repmat(k',1,num_k); %replicate a matrix 1 time for row and 1000 times for column

dis = 1;
tol= 1e-06;
v_guess= zeros(2,num_k);
%%%Consumption function, return utility, iteration%%%
cons_h= A(1)* [kmat .^ alpha]+(1-delta)*kmat - kmat';
returnul_h= cons_h.^(1-sigma)/(1-sigma);
returnul_h(cons_h<0)= -Inf;

cons_l= A(2)* [kmat .^ alpha]+(1-delta)*kmat - kmat';
returnul_l= cons_l.^(1-sigma)/(1-sigma);
returnul_l(cons_l<0)= -Inf;
while dis>tol
        value_mat_h=returnul_h + beta*(tmtx(1,1)* repmat(v_guess(2,:),num_k,1)+tmtx(1,2)*repmat(v_guess(1,:),num_k,1));
        value_mat_l=returnul_l + beta*(tmtx(2,1)* repmat(v_guess(2,:),num_k,1)+tmtx(2,2)*repmat(v_guess(1,:),num_k,1));
        
        [vfnh, pol_indxh]= max(value_mat_h,[],2); %finding largest value along each row (because here we choose k'
        vfnh=vfnh';
        
        [vfnl, pol_indxl]=max(value_mat_l,[],2);
        vfnl=vfnl';
        dis=[max(abs(vfnh- v_guess(2,:)));max(abs(vfnl-v_guess(1,:)))]
        v_guess=[vfnl;vfnh]
end
gh= k(pol_indxh); 
gl=k(pol_indxl);

plot(k,gh,'b',k,gl,'r'),legend('gh','gl')

plot(k,vfnh,'b',k,vfnl,'r'),legend('vfnh','vfnl')









