function Pdot=dos_eq_of_motion_heat_ln_10_12(t,P);
%  DOS_EQ_OF_MOTION : IQT equation of motion 
% --------------------------------------
global number_of_state energy_value number_of_degenerate_energy tau kB beta_reservoir ln_number_of_degenerate_energy
global vq1 iiii rt rs

for I=1:number_of_state
            if(P(I)<10^-26)
                P(I)=10^-26;
            end
end

P=P*(1/sum(P) ); 

sss=0; 
eee=0;
for I=1:number_of_state
    sss=sss-kB*P(I)*(log(P(I))-ln_number_of_degenerate_energy(I));
    eee=eee+P(I)*energy_value(I);
    
    
end
disp(t);
disp(rs);
%vq1(round(t)+1)
for I=1:number_of_state
    pdot(I)=(1/(tau))*P(I)*((kB*-(log(P(I))-ln_number_of_degenerate_energy(I))-sss)-beta_reservoir*rs*(energy_value(I)-eee));
    %pdot(I)=(1/(tau))*P(I)*(-(log(P(I))-ln_number_of_degenerate_energy(I))-sss-beta_reservoir*(energy_value(I)-eee));
end

Pdot=real(pdot');

