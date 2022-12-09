clear all
format longG
close(gcf)
global number_of_state energy_value number_of_degenerate_energy tau Concentration kB
global number_of_degenerate_energy_ln energy_value_ln Concentration_ln inc increment OldP Oldt PP PPP  ln_number_of_degenerate_energy beta_reservoir
global vq1 iiii rt tempSe rs

kB = 8.314;
kB = 1;
amount_solvent=5800;
scale = 1;
q=importdata("DOS.txt");

q(:,2)=q(:,2)-100;
E=(q(:,1));
LNE=log((E(:,1)));
g=exp(q(:,2));
LNg=q(:,2);
ff = realmin;

colorrewl=[ [0,.5,1];[0,.5,1];[0,.5,1];[0,.5,1];[0,.5,1];[0,.5,1];[0,.5,1];...
   [0,.5,1];[0,.5,1]; [178/256,72/256,226/256]; ...
    [178/256,72/256,226/256]; [178/256,72/256,226/256]; [178/256,72/256,226/256];...
   [178/256,72/256,226/256]; [1,0,0]; [1,0,0]; [1,0,0]; [1,0,0]; [1,0,0]; [1,0,0]; [1,0,0]; ...
  [255/256,165/256,0/256];...
     [255/256,165/256,0/256]; [255/256,165/256,0/256]; [255/256,165/256,0/256]; ...
    [255/256,165/256,0/256]; [1,1,0]; [1,1,0];...
    [1,1,0]; [1,1,0]; [1,1,0]; [1,1,0]; [1,1,0];[1,1,0]; [1,1,0]];
plot(E*scale,LNg,'color','black','LineWidth',2')
hold on
scatter(E*scale,LNg,150,'MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor','black','LineWidth',1.5)
ylabel('Density of States, ln(g)'), xlabel('Energy')
set(gca,'FontWeight','bold')
ax = gca;
grid(ax,'on')
set(gca,'GridLineStyle','  -  -')
xtickformat('%.1f')
ytickformat('%.0f')
ax.FontName='Times';
ax.XAxis.FontSize = 16.5;
ax.YAxis.FontSize = 16.5;
ax.LineWidth = 1.5;
%axis([ 0  2.3*10^9  ...
%                0  16000])
C = ["DOS.png"];
str = join(C,"")
exportgraphics(ax ,str,'Resolution',125)


plot(E*scale,LNg,'color','black','LineWidth',2')
hold on
for I=1:length(g(:,1))
scatter(E(I)*scale,LNg(I),150,'MarkerFaceColor',colorrewl(I,:),'MarkerEdgeColor','black','LineWidth',2)
end
ylabel('Density of States, ln(g)'), xlabel('Energy')
set(gca,'FontWeight','bold')
ax = gca;
grid(ax,'on')
set(gca,'GridLineStyle','  -  -')
xtickformat('%.1f')
ytickformat('%.0f')
ax.FontName='Times';
ax.XAxis.FontSize = 16.5;
ax.YAxis.FontSize = 16.5;
ax.LineWidth = 1.5;
%axis([ 0  2.3*10^9  ...
%                0  16000])
C = ["REWLDOS.png"];
str = join(C,"")
exportgraphics(ax ,str,'Resolution',125)


number_of_state = length(q);
tempSe=10;
BetaSe=1/(kB*10);
beta_reservoir = 1/(kB*10);


lnXse=0;

for I=1:number_of_state
    lnXse(I)=LNg(I)-BetaSe*(E(I));
end

pp=max(lnXse);

lnZse=0
for I=1:number_of_state
    lnZse=lnZse+exp(lnXse(I)-pp);
end

lnZse=(pp+log(lnZse));


lnPse=0;
for I=1:number_of_state
    lnPse(I)=LNg(I)-BetaSe*(E(I))-lnZse;
end

Pse=exp(lnPse);

for I=1:number_of_state
    if(Pse(I)<10^-26)
        Pse(I)=10^-26;
    end
end

PseSx = 0;
for U=1:1
    for UU=1:length(Pse(1,:))
        PseSx=PseSx-Pse(U,UU)*(log(Pse(U,UU))-LNg(UU))
    end
end


tempPe=500;
BetaPe=1/(kB*tempPe);

lnXse=0;

for I=1:number_of_state
    lnXpe(I)=LNg(I)-BetaPe*(E(I));
end

pp=max(lnXpe);

lnZpe=0
for I=1:number_of_state
    lnZpe=lnZpe+exp(lnXpe(I)-pp);
end

lnZpe=(pp+log(lnZpe));


lnPpe=0;
for I=1:number_of_state
    lnPpe(I)=LNg(I)-BetaPe*(E(I))-lnZpe;
end

Ppe=exp(lnPpe);

for I=1:number_of_state
    if(Ppe(I)<10^-26)
        Ppe(I)=10^-26;
    end
end

PpeSx = 0;
for U=1:1
    for UU=1:length(Ppe(1,:))
        PpeSx=PpeSx-Ppe(U,UU)*(log(Ppe(U,UU))-LNg(UU))
    end
end


for I=1:number_of_state
    if(Ppe(I)<10^-26)
        Ppe(I)=10^-26;
    end
end



lambda=.0001;

lambda=.00001;
P0=(1-lambda)*Ppe+lambda*Pse;

ExvS=0;
iii=1;
for inc=.001:.02:3000



    BetaSe=1/(kB*inc);
    beta_reservoir = 1/(kB*inc);

    lnXse=0;

    for I=1:number_of_state
        lnXse(I)=LNg(I)-BetaSe*(E(I));
    end

    pp=max(lnXse);

    lnZse=0;
    for I=1:number_of_state
        lnZse=lnZse+exp(lnXse(I)-pp);
    end

    lnZse=(pp+log(lnZse));


    lnPse=0;
    for I=1:number_of_state
        lnPse(I)=LNg(I)-BetaSe*(E(I))-lnZse;
    end

    Pse=exp(lnPse);

    for I=1:number_of_state
        if(Pse(I)<10^-26)
            Pse(I)=10^-26;
        end
    end

    PseSx = 0;
    for U=1:1
        for UU=1:length(Pse(1,:))
            PseSx=PseSx-Pse(U,UU)*(log(Pse(U,UU))-LNg(UU));
        end
    end
    PseEx = 0;
    for U=1:1
        for UU=1:length(Pse(1,:))
            PseEx=PseEx+Pse(U,UU)*E(UU);
        end
    end
    ExvS(iii,1)=PseSx;
    ExvS(iii,2)=PseEx;
    ExvS(iii,3)=inc;
    iii=iii+1;
end
ft=0;
tscale =1;

plot(ExvS(:,1),ExvS(:,2),'color','black','LineWidth',2')
ylabel('<E>'), xlabel('<S>')
set(gca,'FontWeight','bold')
ax = gca;
grid(ax,'on')
set(gca,'GridLineStyle','  -  -')
xtickformat('%.1f')
ytickformat('%.0f')
ax.FontName='Times';
ax.XAxis.FontSize = 16.5;
ax.YAxis.FontSize = 16.5;
ax.LineWidth = 1.5;


%%
for VV=1:1:1

    rs=VV
    mkdir(string(VV))

    for V=0:0


        
       
        output = strcat(string(VV),"/",string(V));
        mkdir(output);

        clear P
        close(gcf)
        close all
        realV=V;
        %V=17
        if(1)


            RR=ExvS(:,1)*scale*kB;
            RR(:,2)=ExvS(:,2)*scale;
            C = [pwd,"/",output,"/","ExvSx",".xlsx"];
            str = join(C,"")
            writematrix(RR,str)


               
                 if(V==0)
                      Temp=5;
                ft=5;

                minI=7;
                maxI=20;
                minII=7;
                maxII=20;
                reducep =30;
                 end

                tspan=[0:.01:2]
                %tspan=[0:50:100000]
            
            
            BetaSe=1/(kB*Temp);
            number_of_state = length(q);
            tempSe=Temp;
            beta_reservoir = 1/(kB*Temp);
            lnXse=0;

            for I=1:number_of_state
                lnXse(I)=LNg(I)-BetaSe*(E(I));
            end

            pp=max(lnXse);

            lnZse=0
            for I=1:number_of_state
                lnZse=lnZse+exp(lnXse(I)-pp);
            end

            lnZse=(pp+log(lnZse));

            lnPse=0;
            for I=1:number_of_state
                lnPse(I)=LNg(I)-BetaSe*(E(I))-lnZse;
            end

            Pse=exp(lnPse);

            for I=1:number_of_state
                if(Pse(I)<10^-26)
                    Pse(I)=10^-26;
                end
            end


            tempPe=ft;
            BetaPe=1/(kB*tempPe);
            for I=1:number_of_state
                lnXpe(I)=LNg(I)-BetaPe*(E(I));
            end
            ppp=max(lnXpe);
            lnZpe=0;
            for I=1:number_of_state
                lnZpe=lnZpe+exp(lnXpe(I)-ppp);
            end
            lnZpe=(ppp+log(lnZpe));
            lnPpe=0;
            for I=1:number_of_state
                lnPpe(I)=LNg(I)-BetaPe*(E(I))-lnZpe;
            end
            Ppe=exp(lnPpe);

           

            tempPe=ft;
            BetaPe=1/(kB*tempPe);
            lnXpe=zeros(1,number_of_state);
            Delta=zeros(1,number_of_state);

            lnXpe=0;
           
            for i=minI:maxI
                Delta(i)=-19000;
            end
            for i=minII:maxII
                Delta(i)=-19000;
            end

%             if(V==0)
%                 Delta=zeros(1,number_of_state);
%                 ft=5;
%             end

            for I=1:number_of_state
                lnXpe(I)=Delta(I)+LNg(I)-BetaPe*(E(I));
            end

            ppp=max(lnXpe);

            lnZpe=0;
            for I=1:number_of_state
                lnZpe=lnZpe+exp(lnXpe(I)-ppp);
            end

            lnZpe=(ppp+log(lnZpe));


            lnPpe=0;

            for I=1:number_of_state
                lnPpe(I)=Delta(I)+LNg(I)-BetaPe*(E(I))-lnZpe;
            end

            Ppe=exp(lnPpe);

            [M,I]=max(Ppe)

            for I=1:number_of_state
                if(Ppe(I)<10^-26)
                    Ppe(I)=10^-26;
                end
            end

            lambda=.001;

            lambda=.001;
            P0=(1-lambda)*Ppe+lambda*Pse;

            for I=1:number_of_state
                if(P0(I)<10^-200)
                    P0(I)=10^-200;
                end
            end

            P0=P0.*1/sum(P0);

            number_of_state = length(E);
            energy_value = E;
            energy_value_ln = LNE;
            number_of_degenerate_energy = g;
            number_of_degenerate_energy_ln = LNg;
            ln_number_of_degenerate_energy = LNg;
            tau=1*tscale
            increment=0;

            lnP0=log(P0);

            ts=0;  %  start time
            tf=5*tau;  %  finish time
            tstep=tf*1*.0001;  %  time step
            %tspan=[0; 30   ]; %GG
            %tspan=[0; 16500 ]; %OR
            %tspan=[0; .02  ]; %OR Smaller
            %tspan=[0; .16   ]; %Sint
            %tspan=[0; 22900  ];  % Better Sint 1214 Guess

            %tspan=[0;.012  ]; %Sint 2



            options=odeset('RelTol',1e-14,'AbsTol',1e-20);
            %options=odeset('RelTol',1e-8,'AbsTol',1e-15);
            [t,P]=ode45('dos_eq_of_motion_heat_ln_10_12',tspan,P0,options);


            for II=1:length(t)
                for I=1:number_of_state
                    if(P(II,I)<10^-26)
                        P(II,I)=10^-26;
                    end
                end
                P(II,:)=P(II,:).*1/sum(P(II,:));
            end


            hold off
            hold off
            % p=plot(t,P,'LineWidth',2);
            % p=plot(t,CompressedP,'LineWidth',2);

            Ex=zeros(1,length(P(:,1)));
            for U=1:length(P(:,1))
                for UU=1:length(P(1,:))
                    Ex(U)=real(Ex(U)+P(U,UU)*E(UU));
                end
            end

            Sx=zeros(1,length(P(:,1)));
            for U=1:length(P(:,1))
                for UU=1:length(P(1,:))
                    Sx(U)=real(Sx(U)-P(U,UU)*(log(P(U,UU))-LNg(UU)));
                end
            end

            plot(ExvS(:,1)*scale*kB,ExvS(:,2)*scale,'Color',[.5 .5 .5],'LineWidth',.5)
            hold on
            plot(Sx*scale*kB,Ex*scale,'color','red','LineWidth',3)
            xlabel('Entropy'), ylabel('Energy')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            axis([ 5  28  ...
                -50  0])
            set(gca,'GridLineStyle','  -  -')
            xtickformat('%.0f')
            ytickformat('%.0f')
            ax.FontName='Times';
            Ex=zeros(1,length(P(:,1)));
            for U=1:length(P(:,1))
                for UU=1:length(P(1,:))
                    Ex(U)=real(Ex(U)+P(U,UU)*E(UU));
                end
            end

            Sx=zeros(1,length(P(:,1)));
            for U=1:length(P(:,1))
                for UU=1:length(P(1,:))
                    Sx(U)=real(Sx(U)-P(U,UU)*(log(P(U,UU))-LNg(UU)));
                end
            end
            plot(ExvS(:,1)*scale*kB,ExvS(:,2)*scale,'Color',[.5 .5 .5],'LineWidth',.5)
            hold on
            plot(Sx*scale*kB,Ex*scale,'color','red','LineWidth',3)
            xlabel('Entropy, J per mol/K'), ylabel('Energy, J/mol')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            set(gca,'GridLineStyle','  -  -')
            xtickformat('%.0f')
            ytickformat('%.0f')
            ax.FontName='Times';

            axis([ 5  28  ...
                -50  0])
            %            hkhk/iouhiohu
            close(gcf)

           
            %         P=P(1:max(EffectiveZeroLoc),:);
            %         t=t(1:max(EffectiveZeroLoc),:);
            %         CompressedP=CompressedP(1:max(EffectiveZeroLoc),:);


            hold off
            hold off
            PPP=turbo(length(P(1,:)));
            
            for I=1:length(P(1,:))
                p=plot(t,P(:,I),'LineWidth',2,'Color',PPP(I,:));
                hold on

                for II=1:reducep:length(P(:,1))
                    
                    mp=max(P(II,:))

                    if(P(II,I)>mp*.1)


                p=scatter(t(II),P(II,I),50,'MarkerEdgeColor','black','MarkerFaceColor',PPP(I,:));
                hold on
                    end
                end
            end
            % p=plot(t,CompressedP,'LineWidth',2)
            %xlabel('Time, s'), ylabel('Eigen Level Probabilities')
            xlabel('time, t*'), ylabel('Energy Level Probabiliy')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            set(gca,'GridLineStyle','  -  -')
            xtickformat('%.1f')
            ytickformat('%.2f')
            ax.FontName='Times';
            ax.XAxis.FontSize = 16.5;
            ax.YAxis.FontSize = 16.5;
            ax.LineWidth = 1.5;


            C = [pwd,"/",output,"/","P",string(realV),"_T_",string(Temp),".png"];
            str = join(C,"")
            exportgraphics(ax ,str,'Resolution',125)

            close(gcf);

                        hold off


            PPP=turbo(length(P(1,:)));

            for L=1:reducep:length(P(:,1))
                %scatter(I,P(L,I));
                %line((1:26),P(L,:))


c = linspace(1,10,length((1:.002:length(q(:,1)))));
qqq=interp1((1:length(q(:,1))),P(L,:),(1:.002:length(q(:,1))));
scatter((1:.002:length(q(:,1))),qqq,5,c,'filled','square');
colormap(gca,"turbo")

                hold on

               % patch([(1:length(q(:,1))) nan],[P(L,:) nan],[P(L,:) nan],[P(L,:) nan],'edgecolor','interp')
               % colormap(jet)
                
            for I=1:length(P(1,:))
                p=scatter(I,P(L,I),50,'MarkerEdgeColor','black','MarkerFaceColor',PPP(I,:),'LineWidth',1.5);
                hold on
            end
            
            end

            % p=plot(t,CompressedP,'LineWidth',2)
            xlabel('Energy Level Index'), ylabel('Energy Level Probability')
           % xlabel('time, s'), ylabel('Energy Level Probabilities')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            set(gca,'GridLineStyle','  -  -')
            xtickformat('%.1f')
            ytickformat('%.2f')
            ax.FontName='Times';
            ax.XAxis.FontSize = 16.5;
            ax.YAxis.FontSize = 16.5;
            ax.LineWidth = 1.5;

            C = [pwd,"/",output,"/","P2",string(realV),"_T_",string(Temp),".png"];
            str = join(C,"")
            exportgraphics(ax ,str,'Resolution',125)

            close(gcf);

            Ex=zeros(1,length(P(:,1)));
            for U=1:length(P(:,1))
                for UU=1:length(P(1,:))
                    Ex(U)=real(Ex(U)+P(U,UU)*E(UU));
                end
            end

            plot(t,Ex*scale,'color','black','LineWidth',2)
            xlabel('time, s'), ylabel('<E>, J/mol')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            set(gca,'GridLineStyle','  -  -')
            xtickformat('%.1f')
            ytickformat('%.1f')
            ax.FontName='Times';
            ax.XAxis.FontSize = 16.5;
            ax.YAxis.FontSize = 16.5;
            ax.LineWidth = 1.5;


            C = [pwd,"/",output,"/","Ex",string(realV),"_T_",string(Temp),".png"];
            str = join(C,"")
            exportgraphics(ax ,str,'Resolution',125)

            plot(t(1:end-1),diff(Ex*scale),'color','black','LineWidth',2)
            xlabel('time, s'), ylabel('d<E>/dt, J/mol*t')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            set(gca,'GridLineStyle','  -  -')
            xtickformat('%.1f')
            ytickformat('%.1f')
            ax.FontName='Times';
            ax.XAxis.FontSize = 16.5;
            ax.YAxis.FontSize = 16.5;
            ax.LineWidth = 1.5;


            C = [pwd,"/",output,"/","ExP",string(realV),"_T_",string(Temp),".png"];
            str = join(C,"")
            exportgraphics(ax ,str,'Resolution',125)

            RR=t;
            RR(:,2)=Ex*scale';
            C = [pwd,"/",output,"/","Ex",string(realV),"_T_",string(Temp),".csv"];
            str = join(C,"")
            writematrix(RR,str)

            RR=t(2:end);
            RR(:,2)=diff(Ex*scale)';
            C = [pwd,"/",output,"/","dEx",string(realV),"_T_",string(Temp),".csv"];
            str = join(C,"")
            writematrix(RR,str)

            Ex=zeros(1,length(P(:,1)));
            Ex2=zeros(1,length(P(:,1)));
            for U=1:length(P(:,1))
                for UU=1:length(P(1,:))
                    Ex(U)=real(Ex(U)+P(U,UU)*E(UU));
                    Ex2(U)=real(Ex2(U)+P(U,UU)*E(UU)^2);
                end
            end

            HC = (Ex2-Ex.^2)/(kB*(Temp*10^-1)^2);
            HC=HC*scale;

            plot(t,HC,'color','black','LineWidth',2)
            xlabel('time, s'), ylabel('Expected Heat Capacity')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            set(gca,'GridLineStyle','  -  -')
            xtickformat('%.1f')
            ytickformat('%.1f')
            ax.FontName='Times';
            ax.XAxis.FontSize = 16.5;
            ax.YAxis.FontSize = 16.5;
            ax.LineWidth = 1.5;


            C = [pwd,"/",output,"/","HC",string(realV),"_T_",string(Temp),".png"];
            str = join(C,"")
            exportgraphics(ax ,str,'Resolution',125)


            Sx=zeros(1,length(P(:,1)));
            for U=1:length(P(:,1))
                for UU=1:length(P(1,:))
                    Sx(U)=real(Sx(U)-P(U,UU)*(log(P(U,UU))-LNg(UU)));
                end
            end

            plot(t,Sx*scale*kB,'color','black','LineWidth',2)
            xlabel('time, s'), ylabel('<S>, J per mol/K')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            set(gca,'GridLineStyle','  -  -')
            xtickformat('%.1f')
            ytickformat('%.1f')
            ax.FontName='Times';
            ax.XAxis.FontSize = 16.5;
            ax.YAxis.FontSize = 16.5;
            ax.LineWidth = 1.5;


            C = [pwd,"/",output,"/","Sx_CycloBrush_",string(realV),"_T_",string(Temp),".png"];
            str = join(C,"")
            exportgraphics(ax ,str,'Resolution',125)

            plot(t,Sx*scale*kB-Sx(1)*scale*kB,'color','black','LineWidth',2)
            xlabel('time, s'), ylabel('Δ<S^A>, J per mol/K')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            set(gca,'GridLineStyle','  -  -')
            xtickformat('%.1f')
            ytickformat('%.1f')
            ax.FontName='Times';
            ax.XAxis.FontSize = 16.5;
            ax.YAxis.FontSize = 16.5;
            ax.LineWidth = 1.5;


            C = [pwd,"/",output,"/","SxA_CycloBrush_",string(realV),"_T_",string(Temp),".png"];
            str = join(C,"")
            exportgraphics(ax ,str,'Resolution',125)

            plot(t(1:end-1),diff(Sx*scale*kB-Sx(1)*scale*kB),'color','black','LineWidth',2)
            xlabel('time, s'), ylabel('d<S>/dt, J per mol/K*t')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            set(gca,'GridLineStyle','  -  -')
            xtickformat('%.1f')
            ytickformat('%.1f')
            ax.FontName='Times';
            ax.XAxis.FontSize = 16.5;
            ax.YAxis.FontSize = 16.5;
            ax.LineWidth = 1.5;


            C = [pwd,"/",output,"/","dSxA_CycloBrush_",string(realV),"_T_",string(Temp),".png"];
            str = join(C,"")
            exportgraphics(ax ,str,'Resolution',125)



            RR=t(1:end);
            RR(:,2)=(Sx*scale*kB-Sx(1)*scale*kB)';
            C = [pwd,"/",output,"/","SxA_CycloBrush_",string(realV),"_T_",string(Temp),".csv"];
            str = join(C,"")
            writematrix(RR,str)

            RR=t(1:end-1);
            RR(:,2)=(diff(Sx*scale*kB-Sx(1)*scale*kB))';
            C = [pwd,"/",output,"/","dSx_CycloBrush_",string(realV),"_T_",string(Temp),".csv"];
            str = join(C,"")
            writematrix(RR,str)

            plot(t(1:end),(Sx*scale*kB-Sx(1)*scale*kB)-(Ex*scale-Ex(1)*scale)/Temp,'color','black','LineWidth',2)
            xlabel('time, s'), ylabel('σ_A, J/K')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            set(gca,'GridLineStyle','  -  -')
            xtickformat('%.1f')
            ytickformat('%.1f')
            ax.FontName='Times';
            ax.XAxis.FontSize = 16.5;
            ax.YAxis.FontSize = 16.5;
            ax.LineWidth = 1.5;


            C = [pwd,"/",output,"/","SxP_CycloBrush_",string(realV),"_T_",string(Temp),".png"];
            str = join(C,"")
            exportgraphics(ax ,str,'Resolution',125)

            RR=t;
            RR(:,2)=((Sx*scale*kB-Sx(1)*scale*kB)-(Ex*scale-Ex(1)*scale)/Temp)';
            C = [pwd,"/",output,"/","SxP_CycloBrush_",string(realV),"_T_",string(Temp),".csv"];
            str = join(C,"")
            writematrix(RR,str)

            plot(t(1:end-1),diff((Sx*scale*kB-Sx(1)*scale*kB)-(Ex*scale-Ex(1)*scale)/Temp),'color','black','LineWidth',2)
            xlabel('time, s'), ylabel('$\dot{\sigma_A}$, J per mol/K*t', 'Interpreter','latex')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            set(gca,'GridLineStyle','  -  -')
            xtickformat('%.1f')
            ytickformat('%.1f')
            ax.FontName='Times';
            ax.XAxis.FontSize = 16.5;
            ax.YAxis.FontSize = 16.5;
            ax.LineWidth = 1.5;


            C = [pwd,"/",output,"/","SigmaP_CycloBrush_",string(realV),"_T_",string(Temp),".png"];
            str = join(C,"")
            exportgraphics(ax ,str,'Resolution',125)

            RR=t(1:end-1);
            RR(:,2)=diff((Sx*scale*kB-Sx(1)*scale*kB)-(Ex*scale-Ex(1)*scale)/Temp)';
            C = [pwd,"/",output,"/","SigmaP_CycloBrush_",string(realV),"_T_",string(Temp),".csv"];
            str = join(C,"")
            writematrix(RR,str)

            plot(ExvS(:,1)*scale*kB,ExvS(:,2)*scale,'black','LineWidth',2)
            hold on
            plot(Sx*scale*kB,Ex*scale,'color','red','LineWidth',2)
            xlabel('Entropy'), ylabel('Energy')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            set(gca,'GridLineStyle','  -  -')
            xtickformat('%0f')
            ytickformat('%.0f')
            ax.FontName='Times';
            axis([ 5  28  ...
                -50  0])
            ax.XAxis.FontSize = 16.5;
            ax.YAxis.FontSize = 16.5;
            ax.LineWidth = 1.5;

            %            knljbkhvghfdg/gyjfthdgx

            hold off
            C = [pwd,"/",output,"/","ExvSx_CycloBrush_",string(realV),"_T_",string(Temp),".png"];
            str = join(C,"")
            exportgraphics(ax ,str,'Resolution',125)

            hold off
            C = [pwd,"/",string(VV),"/","ExvSx_CycloBrush_",string(realV),"_T_",string(Temp),".png"];
            str = join(C,"")
            exportgraphics(ax ,str,'Resolution',125)

            RR=(Sx*scale*kB)';
            RR(:,2)=(Ex*scale)';
            C = [pwd,"/",output,"/","ExvSx_CycloBrush_",string(realV),"_T_",string(Temp),".csv"];
            str = join(C,"")
            writematrix(RR,str)


            plot(ExvS(:,1)*scale*kB,ExvS(:,2)*scale,'Color',[.5 .5 .5],'LineWidth',.5)
            hold on
            plot(Sx*scale*kB,Ex*scale,'color','red','LineWidth',3)
            xlabel('Entropy'), ylabel('Energy')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            set(gca,'GridLineStyle','  -  -')
            xtickformat('%.0f')
            ytickformat('%.0f')
            ax.FontName='Times';

            axis([ 5  28  ...
                -50  0])
            ax.XAxis.FontSize = 16.5;
            ax.YAxis.FontSize = 16.5;
            ax.LineWidth = 1.5;
            %pause(10);

houkbyjvg

            Rgx=zeros(1,length(P(:,1)));
            Tortx=zeros(1,length(P(:,1)));
            Rg2x=zeros(1,length(P(:,1)));
            Tort2x=zeros(1,length(P(:,1)));
            Rgperpx=zeros(1,length(P(:,1)));
            Rgzx=zeros(1,length(P(:,1)));
            SQT=importdata("SysParam.csv");
            Rg=SQT(:,2);
            Tort=SQT(:,3);
            Rg2=SQT(:,7);
            Tort2=SQT(:,8);
            Rgperp=SQT(:,9);
            Rgz=SQT(:,10);

            for U=1:length(P(:,1))
                for UU=1:length(P(1,:))
                    Rgx(U)=Rgx(U)+P(U,UU)*Rg(UU);% Tortuosity
                    Tortx(U)=Tortx(U)+P(U,UU)*Tort(UU);%Radius of Gyration
                    Rg2x(U)=Rg2x(U)+P(U,UU)*Rg2(UU);% Tortuosity
                    Tort2x(U)=Tort2x(U)+P(U,UU)*Tort2(UU);%Radius of Gyration
                    Rgperpx(U)=Rgperpx(U)+P(U,UU)*Rgperp(UU);% Tortuosity
                    Rgzx(U)=Rgzx(U)+P(U,UU)*Rgz(UU);% Tortuosity
                end
            end

            plot(t,(Rg2x),'color','black','LineWidth',2)
            xlabel('time, s'), ylabel('<R_g^2>, Å^2')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            set(gca,'GridLineStyle','  -  -')
            xtickformat('%.1f')
            ytickformat('%.1f')
            ax.FontName='Times';
            ax.XAxis.FontSize = 16.5;
            ax.YAxis.FontSize = 16.5;
            ax.LineWidth = 1.5;


            C = [pwd,"/",output,"/","Rog_CycloBrush_",string(realV),"_T_",string(Temp),".png"];
            str = join(C,"")
            exportgraphics(ax ,str,'Resolution',125)
            hold off
            RR=t;
            RR(:,2)=Rg2x;
            C = [pwd,"/",output,"/","Rog_CycloBrush_",string(realV),"_T_",string(Temp),".csv"];
            str = join(C,"")
            writematrix(RR,str)



            plot(t,(Rg2x),'color','black','LineWidth',2)
            xlabel('time, s'), ylabel('<R_g^2>, Å^2')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            set(gca,'GridLineStyle','  -  -')
            xtickformat('%.1f')
            ytickformat('%.1f')
            ax.FontName='Times';
            ax.XAxis.FontSize = 16.5;
            ax.YAxis.FontSize = 16.5;
            ax.LineWidth = 1.5;



            hold off
            plot(t,(Tort2x),'color','black','LineWidth',2)
            xlabel('time, s'), ylabel('<τ>, Å')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            set(gca,'GridLineStyle','  -  -')
            xtickformat('%.1f')
            ytickformat('%.1f')
            ax.FontName='Times';
            ax.XAxis.FontSize = 16.5;
            ax.YAxis.FontSize = 16.5;
            ax.LineWidth = 1.5;


            C = [pwd,"/",output,"/","Tort_CycloBrush_",string(realV),"_T_",string(Temp),".png"];
            str = join(C,"")
            exportgraphics(ax ,str,'Resolution',125)
            hold off

            RRR=t;
            RRR(:,2)=Tort2x;
            C = [pwd,"/",output,"/","Tort_CycloBrush_",string(realV),"_T_",string(Temp),".csv"];
            str = join(C,"")
            writematrix(RRR,str)

            hold off
            plot(t,(Rgperpx),'color','black','LineWidth',2)
            xlabel('time, s'), ylabel('<R_g_\perp^2>, Å^2')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            set(gca,'GridLineStyle','  -  -')
            xtickformat('%.1f')
            ytickformat('%.1f')
            ax.FontName='Times';
            ax.XAxis.FontSize = 16.5;
            ax.YAxis.FontSize = 16.5;
            ax.LineWidth = 1.5;


            C = [pwd,"/",output,"/","Rogperp_CycloBrush_",string(realV),"_T_",string(Temp),".png"];
            str = join(C,"")
            exportgraphics(ax ,str,'Resolution',125)
            hold off

            RRR=t;
            RRR(:,2)=Rgperpx;
            C = [pwd,"/",output,"/","Rogperp_CycloBrush_",string(realV),"_T_",string(Temp),".csv"];
            str = join(C,"")
            writematrix(RRR,str)

            plot(t(1:end-1),diff(Rgperpx),'color','black','LineWidth',2)
            xlabel('time, s'), ylabel('d<R_g_\perp^2>/dt, Å^2')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            set(gca,'GridLineStyle','  -  -')
            xtickformat('%.1f')
            ytickformat('%.1f')
            ax.FontName='Times';
            ax.XAxis.FontSize = 16.5;
            ax.YAxis.FontSize = 16.5;
            ax.LineWidth = 1.5;


            C = [pwd,"/",output,"/","dRog_CycloBrush_",string(realV),"_T_",string(Temp),".png"];
            str = join(C,"")
            exportgraphics(ax ,str,'Resolution',125)

            RR=t;
            RR(:,2)=Rgx';
            C = [pwd,"/",output,"/","dRog_CycloBrush_",string(realV),"_T_",string(Temp),".csv"];
            str = join(C,"")
            writematrix(RR,str)

            hold off
            plot(t,(Rgzx),'color','black','LineWidth',2)
            xlabel('time, s'), ylabel('<R_g_z^2>, Å^2')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            set(gca,'GridLineStyle','  -  -')
            xtickformat('%.1f')
            ytickformat('%.1f')
            ax.FontName='Times';
            ax.XAxis.FontSize = 16.5;
            ax.YAxis.FontSize = 16.5;
            ax.LineWidth = 1.5;


            C = [pwd,"/",output,"/","Rogz_CycloBrush_",string(realV),"_T_",string(Temp),".png"];
            str = join(C,"")
            exportgraphics(ax ,str,'Resolution',125)
            hold off

            RRR=t;
            RRR(:,2)=Rgzx;
            C = [pwd,"/",output,"/","Rogz_CycloBrush_",string(realV),"_T_",string(Temp),".csv"];
            str = join(C,"")
            writematrix(RRR,str)

            prof=importdata("P_Combined.csv");
            profS=importdata("PS_Combined.csv");
            proff=importdata("P_Combined_Crop.csv");
            proffS=importdata("PS_Combined_Crop.csv");

            reducingindex = 10;

            x=0;

            for xx=1:reducingindex:length(P(:,1))
                x=x+1;
            end

            DPx=zeros(x,length(prof(1,:)));
            DPSx=zeros(x,length(prof(1,:)));
            DPcrx=zeros(x,length(proff(1,:)));
            DPScrx=zeros(x,length(proff(1,:)));
            Heightx=zeros(x,1);
            Widthx=zeros(x,1);

            Abx=zeros(x,1);
            Adx=zeros(x,1);
            Addx=zeros(x,1);
            shortt=zeros(x,1);

            x=1;


            for xx=1:reducingindex:length(P(:,1))

                for UU=1:length(P(1,:))
                    DPx(x,:)=real(DPx(x,:)+prof(UU,:).*P(xx,UU));
                end

                for UU=1:length(P(1,:))
                    DPSx(x,:)=real(DPSx(x,:)+profS(UU,:).*P(xx,UU));
                end

                for UU=1:length(P(1,:))
                    DPcrx(x,:)=real(DPcrx(x,:)+proff(UU,:).*P(xx,UU));
                end

                for UU=1:length(P(1,:))
                    DPScrx(x,:)=real(DPScrx(x,:)+proffS(UU,:).*P(xx,UU));
                end

                ph=1:length(prof(1,:));
                for UU=1:length(P(1,:))
                    plhdrprof=movmean(prof(UU,:),2);
                    Heightx(x)=Heightx(x)+(trapz(ph.*plhdrprof(1,:))/trapz(plhdrprof(1,:))).*P(xx,UU);
                end


                Widthx(x)=0;
                for UU=1:length(P(1,:))
                    plhdrprof=movmean(prof(UU,:),2);
                    ph=gradient(plhdrprof(1,:),[1:80]);
                    max_phi = ph(7);
                    max_DP = plhdrprof(1,:);
                    max_max_DP=max_DP(7);
                    for U=10:length(plhdrprof(1,:))
                        if(ph(U)<max_phi)
                            max_phi = ph(U);
                        end
                        if(max_DP(U)>max_max_DP)
                            max_max_DP = max_DP(U);
                        end
                    end

                    Widthx(x)=Widthx(x)+(max_max_DP/abs(max_phi)).*P(xx,UU);
                end


                Abx(x)=0;
                for UU=1:length(P(1,:))
                    plhdrprofS=movmean(profS(UU,1:end),2);
                    Abx(x)=Abx(x)+((17*17)/(amount_solvent)).*trapz(plhdrprofS(1,1:round(Heightx(x)))).*P(xx,UU);
                end


                moving_index = Heightx(x)+5;
                Adx(x)=0;
                for UU=1:length(P(1,:))
                    plhdrprofS=movmean(profS(UU,1:end),2);
                    gradplhdrprofS=diff(movmean(profS(UU,1:end),2));
                    zeroplhdrprofS=zeros(length(prof(UU,1:end)))+1;



                    lbounds=.982;
                    ubouonds=1.018;
                    for I=1:length(prof(UU,1:end))-7
                        if(plhdrprofS(I+1)*lbounds<plhdrprofS(I) && plhdrprofS(I)<plhdrprofS(I+1)*ubouonds && ...
                                plhdrprofS(I+2)*lbounds<plhdrprofS(I) && plhdrprofS(I)<plhdrprofS(I+2)*ubouonds && ...
                                plhdrprofS(I+3)*lbounds<plhdrprofS(I) && plhdrprofS(I)<plhdrprofS(I+3)*ubouonds )%&& ...
                            % plhdrprofS(I+4)*lbounds<plhdrprofS(I) && plhdrprofS(I)<plhdrprofS(I+4)*ubouonds && ...
                            %  plhdrprofS(I+5)*lbounds<plhdrprofS(I) && plhdrprofS(I)<plhdrprofS(I+5)*ubouonds )
                            zeroplhdrprofS(I) = 0;

                        else
                            zeroplhdrprofS(I) = .6;
                        end

                    end

                    for I=round(Heightx(x))+1:length(prof(UU,1:end))-7
                        if(( zeroplhdrprofS(I)==0 &&...
                                zeroplhdrprofS(I+1)==0 &&...
                                zeroplhdrprofS(I+2)==0 &&...
                                zeroplhdrprofS(I+3)==0 )|| (zeroplhdrprofS(I)==0 &&...
                                zeroplhdrprofS(I-1)==0 &&...
                                zeroplhdrprofS(I-2)==0 &&...
                                zeroplhdrprofS(I-3)==0))%&&...
                            % zeroplhdrprofS(I+4)==0 &&...
                            % zeroplhdrprofS(I+5)==0)

                            moving_index=I;

                            zeroplhdrprofS(I) = -.5;
                            break
                        end

                    end

                    % plot(plhdrprofS(1:end))
                    %  hold on
                    % % %plot( gradplhdrprofS(round(Heightx):round(Heightx)+40))
                    % plot( zeroplhdrprofS)
                    %
                    % % plot(plhdrprofS(round(Heightx):round(Heightx)+40))
                    % % hold on
                    % % %plot( gradplhdrprofS(round(Heightx):round(Heightx)+40))
                    % %  plot( zerogradplhdrprofS(round(Heightx):round(Heightx)+40))



                    %                     prev =  plhdrprofS(1,round(Heightx(x)));
                    %                     moving_index = round(Heightx(x));
                    %                     for U=((round(Heightx(x))):70)
                    %                         if(gradplhdrprofS(1,U)*1000<1 && gradplhdrprofS(1,U+1)*1000<1)
                    %
                    %                         else
                    %                             prev = plhdrprofS(1,U);
                    %                             moving_index=U;
                    %                         end
                    %                     end

                    Adx(x)=Adx(x)+((17*17)/(amount_solvent)).*trapz(plhdrprofS(1,round(Heightx(x)):moving_index)).*P(xx,UU);
                    Addx(x)=Addx(x)+length(round(Heightx(x)):moving_index)*2*((17*17)).*trapz(plhdrprofS(1,round(Heightx(x)):moving_index)).*P(xx,UU);
                end
                x
                shortt(x)=t(xx)+.0000000000001;
                x=x+1;
            end
            close(gcf)
            hold off


            %t2=text(start+width/2,ytick(max(length(ytick))),'//
            % Break_axis_(xx(1:end-10),plhdr(1:end-10),0,15,15)
            % hold on
            % Break_axis_(xx(9:end),plhdr(9:end),2*length(xx)-10,2*length(xx),10)
            % sfsvsfvfs/sfbfsbsfb


            xlabel('Depth, Å'), ylabel('Polymer Brush Density, <ϕ_p>')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            set(gca,'GridLineStyle','  --')
            xtickformat('%0f')
            ytickformat('%.2f')
            ax.FontName='Times';
            ax.XAxis.FontSize = 16.5;
            ax.YAxis.FontSize = 16.5;
            ax.LineWidth = 1.5;

            N=length(DPx(:,1));
            %xx=(1:length(DPSx(1,:)))*2;
            PPP=turbo(N);
            %DPmax=4;
            hold on

            qqqq= movmean(DPx(1,:),2);
            minii=qqqq(1,15)
            mini=1;
            for i=2:length(DPx(:,1))
                qqqq= movmean(DPx(i,:),2);
                qqqq(1,7)
                if(minii>qqqq(1,7))
                    minii= qqqq(1,7);
                    mini=i;
                end
            end

            qqqq= movmean(DPx(1,:),2);
            maxii=max(qqqq);
            maxi=1;
            for i=2:length(DPx(:,1))
                qqqq= movmean(DPx(i,:),2);
                if(maxii>max(qqqq(1,:)))
                    maxii= max(qqqq(1,:));
                    maxi=i;
                end
            end

            maxstart=max(DPx(1,:));
            maxend=max(DPx(end,:));

            for i=1:length(DPx(:,1))
                %         plhdrDP=movmean(DPSx(i,:),2);
                %         if max(plhdrDP)~=0
                %             plot([0,xx(1,DPmax)],[plhdrDP(1,DPmax),plhdrDP(1,DPmax)],':','LineWidth',2,'Color',PPP(i,:))
                %             %plot(xx(1,1:8),plhdrDP(1,1:8))
                %             plot(xx(1,DPmax:69),plhdrDP(1,DPmax:69),'LineWidth',2,'Color',PPP(i,:))
                %             plot([xx(1,69),xx(1,end)],[plhdrDP(1,69),plhdrDP(1,69)],':','LineWidth',2,'Color',PPP(i,:))
                %             i
                %         end

                xx=(1:length(DPx(i,:)))*2.35;
                plhdr = movmean(DPx(i,1:end),2);

                y=plhdr(1:end);
                x=xx(1:end);
                start = 0;
                stop =16;
                width=5;
                %stop =10;
                %width=1;
                y(x>start & x<stop)=[];
                x(x>start & x<stop)=[];
                start2 = 102;
                stop2 = x(end)+1;
                y(x>start2 & x<stop2)=[];
                x(x>start2 & x<stop2)=[];


                % map to new xaxis, leaving a space 'width' wide
                x2=x;
                x2(x2>=stop)=x2(x2>=stop)-(stop-start-width);
                x2(end)=x2(end-1)+x2(end-2)-x2(end-3);
                h=plot(x2,y,'LineWidth',2,'Color',PPP(i,:));
                if(mini==i)
                    ytick=get(gca,'YTick');

                    t11=text(1.2+start+width/2,ytick(1),'▬','fontsize',1,'Color','white');
                    t22=text(1.2+start+width/2,ytick(max(length(ytick))),'▬','fontsize',1,'Color','white');

                    t1=text(start+width/2,ytick(1),'//','fontsize',1);
                    t2=text(start+width/2,ytick(max(length(ytick))),'//','fontsize',1);

                    %                 t11=text(1.2+x2(end),ytick(1),'▬','fontsize',8,'Color','white');
                    %                 t22=text(1.2+x2(end),ytick(max(length(ytick))),'▬','fontsize',8,'Color','white');
                    %
                    %                 t1=text(x2(end),ytick(1),'//','fontsize',18);
                    %                 t2=text(x2(end),ytick(max(length(ytick))),'//','fontsize',18);
                    % For y-axis breaks, use set(t1,'rotation',270);
                    % remap tick marks, and 'erase' them in the gap

                end

                if((maxi==i ))
                    ytick=get(gca,'YTick');

                    %                 t11=text(1.2+start+width/2,ytick(1),'▬','fontsize',8,'Color','white');
                    %                 t22=text(1.2+start+width/2,ytick(max(length(ytick))),'▬','fontsize',8,'Color','white');
                    %
                    %                 t1=text(start+width/2,ytick(1),'//','fontsize',18);
                    %                 t2=text(start+width/2,ytick(max(length(ytick))),'//','fontsize',18);

                    t11=text(1.2+x2(end),ytick(1),'▬','fontsize',1,'Color','white');
                    t22=text(1.2+x2(end),ytick(max(length(ytick))),'▬','fontsize',1,'Color','white');

                    t1=text(x2(end),ytick(1),'//','fontsize',1);
                    t2=text(x2(end),ytick(max(length(ytick))),'//','fontsize',1);
                    % For y-axis breaks, use set(t1,'rotation',270);
                    % remap tick marks, and 'erase' them in the gap

                end

                if(i==1 )
                    xtick=get(gca,'XTick');
                    dtick=xtick(2)-xtick(1);
                    gap=floor(width/dtick);
                    last=max(xtick(xtick<=start));          % last tick mark in LH dataset
                    next=min(xtick(xtick>=(last+dtick*(1+gap))));   % first tick mark within RH dataset
                    offset=size(x2(x2>last&x2<next),2)*(x(2)-x(1));
                    for i=1:sum(xtick>(last+gap))
                        xtick(find(xtick==last)+i+gap)=stop+offset+dtick*(i-1);
                    end
                    xtick(end)=xx(end);
                    for i=1:length(xtick)
                        if xtick(i)>last&xtick(i)<next
                            xticklabel{i}=sprintf('%.1f',[]);
                        else
                            xticklabel{i}=sprintf('%.1f',xtick(i));
                        end
                    end
                    set(gca,'xticklabel',xticklabel);
                end

                if(i==1 )
                    xtick=get(gca,'XTick');
                    dtick=xtick(2)-xtick(1);
                    gap=floor(width/dtick);
                    last=max(xtick(xtick<=start));          % last tick mark in LH dataset
                    next=min(xtick(xtick>=(last+dtick*(1+gap))));   % first tick mark within RH dataset
                    offset=size(x2(x2>last&x2<next),2)*(x(2)-x(1));
                    for i=1:sum(xtick>(last+gap))
                        xtick(find(xtick==last)+i+gap)=stop+offset+dtick*(i-1);
                    end
                    xtick(end)=xx(end);
                    for i=1:length(xtick)
                        if xtick(i)>last&xtick(i)<next
                            xticklabel{i}=sprintf('%.1f',[]);
                        else
                            xticklabel{i}=sprintf('%.1f',xtick(i));
                        end
                    end
                    set(gca,'xticklabel',xticklabel);
                end

            end

            for i=1:length(DPx(:,1))
                %         plhdrDP=movmean(DPSx(i,:),2);
                %         if max(plhdrDP)~=0
                %             plot([0,xx(1,DPmax)],[plhdrDP(1,DPmax),plhdrDP(1,DPmax)],':','LineWidth',2,'Color',PPP(i,:))
                %             %plot(xx(1,1:8),plhdrDP(1,1:8))
                %             plot(xx(1,DPmax:69),plhdrDP(1,DPmax:69),'LineWidth',2,'Color',PPP(i,:))
                %             plot([xx(1,69),xx(1,end)],[plhdrDP(1,69),plhdrDP(1,69)],':','LineWidth',2,'Color',PPP(i,:))
                %             i
                %         end

                xx=(1:length(DPx(i,:)))*2.35;
                plhdr = movmean(DPx(i,1:end),2);

                y=plhdr(1:end);
                x=xx(1:end);
                start = 0;
                stop =16;
                width=5;
                %stop =10;
                %width=1;
                y(x>start & x<stop)=[];
                x(x>start & x<stop)=[];
                start2 = 102;
                stop2 = x(end)+1;
                y(x>start2 & x<stop2)=[];
                x(x>start2 & x<stop2)=[];


                % map to new xaxis, leaving a space 'width' wide
                x2=x;
                x2(x2>=stop)=x2(x2>=stop)-(stop-start-width);
                x2(end)=x2(end-1)+x2(end-2)-x2(end-3);
                h=plot(x2,y,'LineWidth',2,'Color',PPP(i,:));
                if(mini==i)
                    ytick=get(gca,'YTick');

                    t11=text(1.2+start+width/2,ytick(1),'▬','fontsize',8,'Color','white');
                    t22=text(1.2+start+width/2,ytick(max(length(ytick))),'▬','fontsize',8,'Color','white');

                    t1=text(start+width/2,ytick(1),'//','fontsize',18);
                    t2=text(start+width/2,ytick(max(length(ytick))),'//','fontsize',18);

                    %                 t11=text(1.2+x2(end),ytick(1),'▬','fontsize',8,'Color','white');
                    %                 t22=text(1.2+x2(end),ytick(max(length(ytick))),'▬','fontsize',8,'Color','white');
                    %
                    %                 t1=text(x2(end),ytick(1),'//','fontsize',18);
                    %                 t2=text(x2(end),ytick(max(length(ytick))),'//','fontsize',18);
                    % For y-axis breaks, use set(t1,'rotation',270);
                    % remap tick marks, and 'erase' them in the gap

                end

                if((maxi==i ))
                    ytick=get(gca,'YTick');

                    %                 t11=text(1.2+start+width/2,ytick(1),'▬','fontsize',8,'Color','white');
                    %                 t22=text(1.2+start+width/2,ytick(max(length(ytick))),'▬','fontsize',8,'Color','white');
                    %
                    %                 t1=text(start+width/2,ytick(1),'//','fontsize',18);
                    %                 t2=text(start+width/2,ytick(max(length(ytick))),'//','fontsize',18);

                    t11=text(1.2+x2(end),ytick(1),'▬','fontsize',8,'Color','white');
                    t22=text(1.2+x2(end),ytick(max(length(ytick))),'▬','fontsize',8,'Color','white');

                    t1=text(x2(end),ytick(1),'//','fontsize',18);
                    t2=text(x2(end),ytick(max(length(ytick))),'//','fontsize',18);
                    % For y-axis breaks, use set(t1,'rotation',270);
                    % remap tick marks, and 'erase' them in the gap

                end

                if(i==1 )
                    xtick=get(gca,'XTick');
                    dtick=xtick(2)-xtick(1);
                    gap=floor(width/dtick);
                    last=max(xtick(xtick<=start));          % last tick mark in LH dataset
                    next=min(xtick(xtick>=(last+dtick*(1+gap))));   % first tick mark within RH dataset
                    offset=size(x2(x2>last&x2<next),2)*(x(2)-x(1));
                    for i=1:sum(xtick>(last+gap))
                        xtick(find(xtick==last)+i+gap)=stop+offset+dtick*(i-1);
                    end
                    xtick(end)=xx(end);
                    for i=1:length(xtick)
                        if xtick(i)>last&xtick(i)<next
                            xticklabel{i}=sprintf('%.1f',[]);
                        else
                            xticklabel{i}=sprintf('%.1f',xtick(i));
                        end
                    end
                    set(gca,'xticklabel',xticklabel);
                end

                if(i==1 )
                    xtick=get(gca,'XTick');
                    dtick=xtick(2)-xtick(1);
                    gap=floor(width/dtick);
                    last=max(xtick(xtick<=start));          % last tick mark in LH dataset
                    next=min(xtick(xtick>=(last+dtick*(1+gap))));   % first tick mark within RH dataset
                    offset=size(x2(x2>last&x2<next),2)*(x(2)-x(1));
                    for i=1:sum(xtick>(last+gap))
                        xtick(find(xtick==last)+i+gap)=stop+offset+dtick*(i-1);
                    end
                    xtick(end)=xx(end);
                    for i=1:length(xtick)
                        if xtick(i)>last&xtick(i)<next
                            xticklabel{i}=sprintf('%.1f',[]);
                        else
                            xticklabel{i}=sprintf('%.1f',xtick(i));
                        end
                    end
                    set(gca,'xticklabel',xticklabel);
                end

            end
            hold off


            C = [pwd,"/",output,"/","Density_Poly_CycloBrush_",string(realV),"_T_",string(Temp),".png"];
            str = join(C,"")
            exportgraphics(ax ,str,'Resolution',125)

            hold off


            hold off
            close(gcf)
            xlabel('Depth, Å'), ylabel('Solvent Density , <ϕ_s>')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            set(gca,'GridLineStyle','  -  -')
            xtickformat('%0f')
            ytickformat('%.2f')
            ax.FontName='Times';
            ax.XAxis.FontSize = 16.5;
            ax.YAxis.FontSize = 16.5;
            ax.LineWidth = 1.5;

            hold on
            %     plot(movmean(DPSx(1,:),2))
            %     for i=1:length(DPx(:,1))
            %         plot(movmean(DPSx(i,:),2))
            %     end
            %
            N=length(DPSx(:,1));
            %xx=(1:length(DPSx(1,:)))*2;
            PPP=turbo(N);
            %DPmax=4;
            hold on

            qqqq= movmean(DPSx(1,:),2);
            minii=qqqq(1,7)
            mini=1;
            for i=2:length(DPSx(:,1))
                qqqq= movmean(DPSx(i,:),2);
                qqqq(1,7)
                if(minii>qqqq(1,7))
                    minii= qqqq(1,7);
                    mini=i;
                end
            end

            qqqq= movmean(DPSx(1,:),2);
            maxii=max(qqqq);
            maxi=1;
            for i=2:length(DPSx(:,1))
                qqqq= movmean(DPSx(i,:),2);
                if(maxii>max(qqqq(1,:)))
                    maxii= max(qqqq(1,:));
                    maxi=i;
                end
            end

            maxstart=max(DPSx(1,:));
            maxend=max(DPSx(end,:));

            for i=1:length(DPSx(:,1))
                %         plhdrDP=movmean(DPSx(i,:),2);
                %         if max(plhdrDP)~=0
                %             plot([0,xx(1,DPmax)],[plhdrDP(1,DPmax),plhdrDP(1,DPmax)],':','LineWidth',2,'Color',PPP(i,:))
                %             %plot(xx(1,1:8),plhdrDP(1,1:8))
                %             plot(xx(1,DPmax:69),plhdrDP(1,DPmax:69),'LineWidth',2,'Color',PPP(i,:))
                %             plot([xx(1,69),xx(1,end)],[plhdrDP(1,69),plhdrDP(1,69)],':','LineWidth',2,'Color',PPP(i,:))
                %             i
                %         end

                xx=(1:length(DPSx(i,:)))*2.35;
                plhdr = movmean(DPSx(i,1:end),2);

                y=plhdr(1:end);
                x=xx(1:end);
                start = 0;
                stop =8;
                width=6;
                %stop =10;
                %width=1;
                y(x>start & x<stop)=[];
                x(x>start & x<stop)=[];
                start2 = 102;
                stop2 = x(end)+1;
                y(x>start2 & x<stop2)=[];
                x(x>start2 & x<stop2)=[];


                % map to new xaxis, leaving a space 'width' wide
                x2=x;
                x2(x2>=stop)=x2(x2>=stop)-(stop-start-width);
                x2(end)=x2(end-1)+x2(end-2)-x2(end-3);
                h=plot(x2,y,'LineWidth',2,'Color',PPP(i,:));
                if(mini==i)
                    ytick=get(gca,'YTick');

                    t11=text(1.2+start+width/2,ytick(1),'▬','fontsize',8,'Color','white');
                    t22=text(1.2+start+width/2,ytick(max(length(ytick))),'▬','fontsize',8,'Color','white');

                    t1=text(start+width/2,ytick(1),'//','fontsize',18);
                    t2=text(start+width/2,ytick(max(length(ytick))),'//','fontsize',18);

                    %                 t11=text(1.2+x2(end),ytick(1),'▬','fontsize',8,'Color','white');
                    %                 t22=text(1.2+x2(end),ytick(max(length(ytick))),'▬','fontsize',8,'Color','white');
                    %
                    %                 t1=text(x2(end),ytick(1),'//','fontsize',18);
                    %                 t2=text(x2(end),ytick(max(length(ytick))),'//','fontsize',18);
                    % For y-axis breaks, use set(t1,'rotation',270);
                    % remap tick marks, and 'erase' them in the gap

                end

                if((maxi==i ))
                    ytick=get(gca,'YTick');

                    %                 t11=text(1.2+start+width/2,ytick(1),'▬','fontsize',8,'Color','white');
                    %                 t22=text(1.2+start+width/2,ytick(max(length(ytick))),'▬','fontsize',8,'Color','white');
                    %
                    %                 t1=text(start+width/2,ytick(1),'//','fontsize',18);
                    %                 t2=text(start+width/2,ytick(max(length(ytick))),'//','fontsize',18);

                    t11=text(1.2+x2(end),ytick(1),'▬','fontsize',8,'Color','white');
                    t22=text(1.2+x2(end),ytick(max(length(ytick))),'▬','fontsize',8,'Color','white');

                    t1=text(x2(end),ytick(1),'//','fontsize',18);
                    t2=text(x2(end),ytick(max(length(ytick))),'//','fontsize',18);
                    % For y-axis breaks, use set(t1,'rotation',270);
                    % remap tick marks, and 'erase' them in the gap

                end

                if(i==1 )
                    xtick=get(gca,'XTick');
                    dtick=xtick(2)-xtick(1);
                    gap=floor(width/dtick);
                    last=max(xtick(xtick<=start));          % last tick mark in LH dataset
                    next=min(xtick(xtick>=(last+dtick*(1+gap))));   % first tick mark within RH dataset
                    offset=size(x2(x2>last&x2<next),2)*(x(2)-x(1));
                    for i=1:sum(xtick>(last+gap))
                        xtick(find(xtick==last)+i+gap)=stop+offset+dtick*(i-1);
                    end
                    xtick(end)=xx(end);
                    for i=1:length(xtick)
                        if xtick(i)>last&xtick(i)<next
                            xticklabel{i}=sprintf('%.1f',[]);
                        else
                            xticklabel{i}=sprintf('%.1f',xtick(i));
                        end
                    end
                    set(gca,'xticklabel',xticklabel);
                end

            end
            hold off

            C = [pwd,"/",output,"/","Density_Solvent_CycloBrush_",string(realV),"_T_",string(Temp),".png"];
            str = join(C,"")
            exportgraphics(ax ,str,'Resolution',125)

            hold off
            %         plhdr = DPcrx;
            %         pdly = DPx;
            %         for i=1:length(DPcrx(:,1))
            %             pdly(i,:)=movmean(DPx(i,:),2);
            %             plhldr(i,:)=pdly(i,6:75);
            %         end
            %         plot(movmean(plhldr(1,:),1))
            %         xlabel('z axis-depth'), ylabel('Brush Density')
            %         set(gca,'FontWeight','bold')
            %         ax = gca;
            %         grid(ax,'on')
            %         set(gca,'GridLineStyle','  -  -')
            %         xtickformat('%.2f')
            %         ytickformat('%.2f')
            %         ax.FontName='Times';
            %         ax.XAxis.FontSize = 16.5;
            %         ax.YAxis.FontSize = 16.5;
            %         ax.LineWidth = 1.5;
            %
            %         hold on
            %         for i=1:length(plhldr(:,1))
            %             plot(movmean(plhldr(i,:),1))
            %         end
            %         hold off
            %
            %         C = [pwd,"/",output,"/","Density_Poly_CycloBrushCrop_",string(realV),"_T_",string(Temp),".png"];
            %         str = join(C,"")
            %         exportgraphics(ax ,str,'Resolution',125)
            %
            %         hold off
            %
            %         hold off
            %
            %         xlabel('z axis-depth'), ylabel('Solvent Density')
            %         set(gca,'FontWeight','bold')
            %         ax = gca;
            %         grid(ax,'on')
            %         set(gca,'GridLineStyle','  -  -')
            %         xtickformat('%.2f')
            %         ytickformat('%.2f')
            %         ax.FontName='Times';
            %         ax.XAxis.FontSize = 16.5;
            %         ax.YAxis.FontSize = 16.5;
            %         ax.LineWidth = 1.5;
            %
            %         hold on
            %         plhdr = DPScrx;
            %         pdly = DPSx;
            %         for i=1:length(DPScrx(:,1))
            %             pdly(i,:)=movmean(DPSx(i,:),2);
            %             plhldr(i,:)=pdly(i,6:75);
            %         end
            %         plot(movmean(plhldr(1,:),1))
            %         for i=1:length(DPcrx(:,1))
            %             plot(movmean(plhldr(1,:),1))
            %         end
            %         hold off
            %
            %         C = [pwd,"/",output,"/","Density_Solvent_CycloBrushCrop_",string(realV),"_T_",string(Temp),".png"];
            %         str = join(C,"")
            %         exportgraphics(ax ,str,'Resolution',125)
            %
            %         hold off

            plot(shortt,(Heightx)*2.35,'color','black','LineWidth',2)
            xlabel('time, s'), ylabel('Brush Height, <h>, Å')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            set(gca,'GridLineStyle','  -  -')
            xtickformat('%.2f')
            ytickformat('%.2f')
            ax.FontName='Times';
            ax.XAxis.FontSize = 16.5;
            ax.YAxis.FontSize = 16.5;
            ax.LineWidth = 1.5;

            ax.FontName='Times';
            ax.XAxis.FontSize = 16;
            ax.YAxis.FontSize = 16;


            hold off
            plot(shortt,(Heightx)*2.35,'color','black','LineWidth',2)
            xlabel('time, s'), ylabel('Expected Brush Height, <h>, Å')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            set(gca,'GridLineStyle','  -  -')
            xtickformat('%.2f')
            ytickformat('%.2f')
            ax.FontName='Times';
            ax.XAxis.FontSize = 16.5;
            ax.YAxis.FontSize = 16.5;
            ax.LineWidth = 1.5;


            C = [pwd,"/",output,"/","Height_CycloBrush_",string(realV),"_T_",string(Temp),".png"];
            str = join(C,"")
            exportgraphics(ax ,str,'Resolution',125)
            hold off

            RRR=shortt;
            RRR(:,2)=Heightx;
            C = [pwd,"/",output,"/","Height_CycloBrush_",string(realV),"_T_",string(Temp),".csv"];
            str = join(C,"")
            writematrix(RRR,str)

            hold off
            plot(shortt,(Widthx),'color','black','LineWidth',2)
            xlabel('time, s'), ylabel('Expected Interfacial Width, <W>,Å*dϕ_p/dÅ')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            xtickformat('%.2f')
            ytickformat('%.2f')
            ax.FontName='Times';
            ax.XAxis.FontSize = 16;
            ax.YAxis.FontSize = 16;

            C = [pwd,"/",output,"/","Width_CycloBrush_",string(realV),"_T_",string(Temp),".png"];
            str = join(C,"")
            exportgraphics(ax ,str,'Resolution',125)
            hold off

            RRR=shortt;
            RRR(:,2)=Widthx';
            C = [pwd,"/",output,"/","Width_CycloBrush_",string(realV),"_T_",string(Temp),".csv"];
            str = join(C,"")
            writematrix(RRR,str)

            hold off
            plot(shortt,(Abx)*100,'color','black','LineWidth',2)
            xlabel('time, s'), ylabel('Expected Absorbance, <%_a_b>')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            xtickformat('%.2f')
            ytickformat('%.2f')
            ax.FontName='Times';
            ax.XAxis.FontSize = 16;
            ax.YAxis.FontSize = 16;

            C = [pwd,"/",output,"/","Absorbance_CycloBrush_",string(realV),"_T_",string(Temp),".png"];
            str = join(C,"")
            exportgraphics(ax ,str,'Resolution',125)
            hold off

            RRR=shortt;
            RRR(:,2)=Abx*100;
            C = [pwd,"/",output,"/","Absorbance_CycloBrush_",string(realV),"_T_",string(Temp),".csv"];
            str = join(C,"")
            writematrix(RRR,str)

            hold off
            plot(shortt,(Adx)*100,'color','black','LineWidth',2)
            xlabel('time, s'), ylabel('Expected Adsorbance, <%_a_d>')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            xtickformat('%.2f')
            ytickformat('%.1f')
            ax.FontName='Times';
            ax.XAxis.FontSize = 16;
            ax.YAxis.FontSize = 16;

            C = [pwd,"/",output,"/","Adsorbance_CycloBrush_",string(realV),"_T_",string(Temp),".png"];
            str = join(C,"")
            exportgraphics(ax ,str,'Resolution',125)
            hold off

            RRR=shortt;
            RRR(:,2)=Adx*100;
            C = [pwd,"/",output,"/","Adsorbance_CycloBrush_",string(realV),"_T_",string(Temp),".csv"];
            str = join(C,"")
            writematrix(RRR,str)

            hold off
            plot(shortt,(Addx),'color','black','LineWidth',2)
            xlabel('time, s'), ylabel('Expected Adsorbance, <%_a_d> Å^2')
            set(gca,'FontWeight','bold')
            ax = gca;
            grid(ax,'on')
            xtickformat('%.2f')
            ytickformat('%.1f')
            ax.FontName='Times';
            ax.XAxis.FontSize = 16;
            ax.YAxis.FontSize = 16;

            C = [pwd,"/",output,"/","Also_Adsorbance_CycloBrush_",string(realV),"_T_",string(Temp),".png"];
            str = join(C,"")
            exportgraphics(ax ,str,'Resolution',125)
            hold off

            RRR=shortt;
            RRR(:,2)=Addx*100;
            C = [pwd,"/",output,"/","Also_Adsorbance_CycloBrush_",string(realV),"_T_",string(Temp),".csv"];
            str = join(C,"")
            writematrix(RRR,str)
            hold off

            if(SubValue(1)~=0)

                for L=1:length(SubValue)
                    NewP=P;
                    Newt=t;
                    NewP=P(1:SubValue(L),:);
                    Newt=t(1:SubValue(L),:);
                    NewCompressedP=CompressedP(1:SubValue(L),:);


                    hold off
                    hold off
                    % p=plot(t,P,'LineWidth',2);
                    %p=plot(Newt,NewCompressedP,'LineWidth',2)
                    xlabel('Dimensionless Time'), ylabel('Eigen Level Probabilities')
                    xlabel('time, s'), ylabel('Energy Level Probabilities')
                    set(gca,'FontWeight','bold')
                    ax = gca;
                    grid(ax,'on')
                    set(gca,'GridLineStyle','  -  -')
                    xtickformat('%.1f')
                    ytickformat('%.1f')
                    ax.FontName='Times';
                    ax.XAxis.FontSize = 16.5;
                    ax.YAxis.FontSize = 16.5;
                    ax.LineWidth = 1.5;

                    %                 C = [pwd,"/",output,"/","Prob_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".png"];
                    %                 str = join(C,"")
                    %                 exportgraphics(ax ,str,'Resolution',125)
                    %
                    %                 C = [pwd,"/",output,"/","CProb_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".csv"];
                    %                 str = join(C,"")
                    %                 writematrix(CompressedP,str)
                    %C = ["CTime_CycloBrush_",string(realV),"_T_",string(Temp),".csv"];
                    %str = join(C,"")
                    %writematrix(CompressedP,str)

                    close(gcf)
                    close all

                    plot(Newt,Ex(1:SubValue(L))*scale,'color','black','LineWidth',2)
                    xlabel('time, s'), ylabel('<E>, J/mol')
                    set(gca,'FontWeight','bold')
                    ax = gca;
                    grid(ax,'on')
                    set(gca,'GridLineStyle','  -  -')
                    xtickformat('%.1f')
                    ytickformat('%.1f')
                    ax.FontName='Times';
                    ax.XAxis.FontSize = 16.5;
                    ax.YAxis.FontSize = 16.5;
                    ax.LineWidth = 1.5;


                    C = [pwd,"/",output,"/","Ex_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".png"];
                    str = join(C,"")
                    exportgraphics(ax ,str,'Resolution',125)

                    plot(Newt(1:end-1),diff(Ex(1:SubValue(L))*scale),'color','black','LineWidth',2)
                    xlabel('time, s'), ylabel('d<E>/dt*, J/mol*t*')
                    set(gca,'FontWeight','bold')
                    ax = gca;
                    grid(ax,'on')
                    set(gca,'GridLineStyle','  -  -')
                    xtickformat('%.1f')
                    ytickformat('%.1f')
                    ax.FontName='Times';
                    ax.XAxis.FontSize = 16.5;
                    ax.YAxis.FontSize = 16.5;
                    ax.LineWidth = 1.5;


                    C = [pwd,"/",output,"/","ExP_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".png"];
                    str = join(C,"")
                    exportgraphics(ax ,str,'Resolution',125)

                    RR=t(1:SubValue(L));
                    RR(:,2)=(Ex(1:SubValue(L))*scale)';
                    C = [pwd,"/",output,"/","Ex_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".csv"];
                    str = join(C,"")
                    writematrix(RR,str)


                    RR=t(2:SubValue(L));
                    RR(:,2)=diff(Ex(1:SubValue(L))*scale)';
                    C = [pwd,"/",output,"/","dEx_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".csv"];
                    str = join(C,"")
                    writematrix(RR,str)

                    HC = (Ex2(1:SubValue(L))-Ex(1:SubValue(L)).^2)/(kB*(Temp)^2);
                    HC=HC*scale;

                    plot(Newt,HC,'color','black','LineWidth',2)
                    xlabel('time, s'), ylabel('Expected Heat Capacity, J/mol')
                    set(gca,'FontWeight','bold')
                    ax = gca;
                    grid(ax,'on')
                    set(gca,'GridLineStyle','  -  -')
                    xtickformat('%.1f')
                    ytickformat('%.1f')
                    ax.FontName='Times';
                    ax.XAxis.FontSize = 16.5;
                    ax.YAxis.FontSize = 16.5;
                    ax.LineWidth = 1.5;


                    C = [pwd,"/",output,"/","HC_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".png"];
                    str = join(C,"")
                    exportgraphics(ax ,str,'Resolution',125)



                    plot(Newt,Sx(1:SubValue(L))*scale*kB,'color','black','LineWidth',2)
                    xlabel('time, s'), ylabel('<S>, J per mol/K')
                    set(gca,'FontWeight','bold')
                    ax = gca;
                    grid(ax,'on')
                    set(gca,'GridLineStyle','  -  -')
                    xtickformat('%.1f')
                    ytickformat('%.1f')
                    ax.FontName='Times';
                    ax.XAxis.FontSize = 16.5;
                    ax.YAxis.FontSize = 16.5;
                    ax.LineWidth = 1.5;


                    C = [pwd,"/",output,"/","Sx_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".png"];
                    str = join(C,"")
                    exportgraphics(ax ,str,'Resolution',125)

                    plot(Newt,Sx(1:SubValue(L))*scale*kB-Sx(1)*scale*kB,'color','black','LineWidth',2)
                    xlabel('time, s'), ylabel('<S^A>, J per mol/K')
                    set(gca,'FontWeight','bold')
                    ax = gca;
                    grid(ax,'on')
                    set(gca,'GridLineStyle','  -  -')
                    xtickformat('%.1f')
                    ytickformat('%.1f')
                    ax.FontName='Times';
                    ax.XAxis.FontSize = 16.5;
                    ax.YAxis.FontSize = 16.5;
                    ax.LineWidth = 1.5;


                    C = [pwd,"/",output,"/","SxA_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".png"];
                    str = join(C,"")
                    exportgraphics(ax ,str,'Resolution',125)

                    RR=Newt;
                    RR(:,2)=(Sx(1:SubValue(L))*scale*kB-Sx(1)*scale*kB)';
                    C = [pwd,"/",output,"/","SxA_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".csv"];
                    str = join(C,"")
                    writematrix(RR,str)

                    plot(Newt(1:end-1),diff(Sx(1:SubValue(L))*scale*kB),'color','black','LineWidth',2)
                    xlabel('time, s'), ylabel('d<S>/dt*, J per mol/K*t*')
                    set(gca,'FontWeight','bold')
                    ax = gca;
                    grid(ax,'on')
                    set(gca,'GridLineStyle','  -  -')
                    xtickformat('%.1f')
                    ytickformat('%.1f')
                    ax.FontName='Times';
                    ax.XAxis.FontSize = 16.5;
                    ax.YAxis.FontSize = 16.5;
                    ax.LineWidth = 1.5;


                    C = [pwd,"/",output,"/","dSx_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".png"];
                    str = join(C,"")
                    exportgraphics(ax ,str,'Resolution',125)

                    RR=Newt(1:end-1);
                    RR(:,2)=(diff(Sx(1:SubValue(L))*scale*kB))';
                    C = [pwd,"/",output,"/","dSx_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".csv"];
                    str = join(C,"")
                    writematrix(RR,str)

                    plot(Newt(1:end),(Sx(1:SubValue(L))*scale*kB-Sx(1)*scale*kB)-(Ex(1:SubValue(L))*scale-Ex(1)*scale)/Temp,'color','black','LineWidth',2)
                    xlabel('time, s'), ylabel('Δ<S_A>, J/K')
                    set(gca,'FontWeight','bold')
                    ax = gca;
                    grid(ax,'on')
                    set(gca,'GridLineStyle','  -  -')
                    xtickformat('%.1f')
                    ytickformat('%.1f')
                    ax.FontName='Times';
                    ax.XAxis.FontSize = 16.5;
                    ax.YAxis.FontSize = 16.5;
                    ax.LineWidth = 1.5;


                    C = [pwd,"/",output,"/","SxP_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".png"];
                    str = join(C,"")
                    exportgraphics(ax ,str,'Resolution',125)

                    RR=Newt(1:end);
                    RR(:,2)=((Sx(1:SubValue(L))*scale*kB-Sx(1)*scale*kB)-(Ex(1:SubValue(L))*scale-Ex(1)*scale)/Temp)';
                    C = [pwd,"/",output,"/","SxP_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".csv"];
                    str = join(C,"")
                    writematrix(RR,str)

                    plot(Newt(1:end-1),diff((Sx(1:SubValue(L))*scale*kB-Sx(1)*scale*kB)-(Ex(1:SubValue(L))*scale-Ex(1)*scale)/Temp),'color','black','LineWidth',2)
                    xlabel('time, s'), ylabel('$\dot{\sigma_A}$, J per mol/K*t', 'Interpreter','latex')
                    set(gca,'FontWeight','bold')
                    ax = gca;
                    grid(ax,'on')
                    set(gca,'GridLineStyle','  -  -')
                    xtickformat('%.1f')
                    ytickformat('%.1f')
                    ax.FontName='Times';
                    ax.XAxis.FontSize = 16.5;
                    ax.YAxis.FontSize = 16.5;
                    ax.LineWidth = 1.5;


                    C = [pwd,"/",output,"/","SigmaP_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".png"];
                    str = join(C,"")
                    exportgraphics(ax ,str,'Resolution',125)

                    RR=Newt(1:end-1);
                    RR(:,2)=diff((Sx(1:SubValue(L))*scale*kB-Sx(1)*scale*kB)-(Ex(1:SubValue(L))*scale-Ex(1)*scale)/Temp)';
                    C = [pwd,"/",output,"/","SigmaP_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".csv"];
                    str = join(C,"")
                    writematrix(RR,str)

                    plot(ExvS(1:end,1)*scale*kB,ExvS(1:end,2)*scale,'black','LineWidth',2)
                    hold on
                    plot(Sx(1:SubValue(L))*scale*kB,Ex(1:SubValue(L))*scale,'color','red','LineWidth',2)
                    xlabel('Entropy, J per mol/K'), ylabel('Energy, J/mol')
                    set(gca,'FontWeight','bold')
                    ax = gca;
                    grid(ax,'on')
                    set(gca,'GridLineStyle','  -  -')
                    xtickformat('%.1f')
                    ytickformat('%.1f')
                    ax.FontName='Times';
                    axis([ (min(Sx*scale)-.6*10^6)*kB  (max(Sx*scale)+.2*10^6)*kB  ...
                        (min(Ex*scale)-.8*10^7)*1  max(Ex*scale)+.5*10^7])
                    ax.XAxis.FontSize = 16.5;
                    ax.YAxis.FontSize = 16.5;
                    ax.LineWidth = 1.5;

                    hold off
                    C = [pwd,"/",output,"/","ExvSx_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".png"];
                    str = join(C,"")
                    exportgraphics(ax ,str,'Resolution',125)

                    RR=(Sx(1:SubValue(L))*scale*kB)';
                    RR(:,2)=(Ex(1:SubValue(L))*scale)';
                    C = [pwd,"/",output,"/","ExvSx_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".csv"];
                    str = join(C,"")
                    writematrix(RR,str)


                    close(gcf)


                    plot(ExvS(1:end,1)*scale*kB+.02*10^6*kB,ExvS(1:end,2)*scale-.05*10^7,'black','LineWidth',4)
                    hold on
                    plot(ExvS(1:end,1)*scale*kB,ExvS(1:end,2)*scale,'Color',[.5 .5 .5],'LineWidth',.5)
                    hold on
                    plot(Sx(1:SubValue(L))*scale*kB,Ex(1:SubValue(L))*scale,'color','red','LineWidth',3)
                    xlabel('Entropy, J per mol/K'), ylabel('Energy, J/mol')
                    set(gca,'FontWeight','bold')
                    ax = gca;
                    grid(ax,'on')
                    set(gca,'GridLineStyle','  -  -')
                    xtickformat('%.1f')
                    ytickformat('%.1f')
                    ax.FontName='Times';

                    axis([ (min(Sx*scale)-.6*10^6 )*kB (max(Sx*scale)+.2*10^6 )*kB ...
                        min(Ex*scale)-.8*10^7  max(Ex*scale)+.5*10^7])

                    ax.XAxis.FontSize = 16.5;
                    ax.YAxis.FontSize = 16.5;


                    hold off
                    C = [pwd,"/",output,"/","ExvSx_ext_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".png"];
                    str = join(C,"")
                    exportgraphics(ax ,str,'Resolution',125)

                    %pause(10);

                    close(gcf)
                    plot(Newt,(Rg2x(1:SubValue(L))),'color','black','LineWidth',2)
                    xlabel('time, s'), ylabel('<R_g^2>, Å^2')
                    set(gca,'FontWeight','bold')
                    ax = gca;
                    grid(ax,'on')
                    set(gca,'GridLineStyle','  -  -')
                    xtickformat('%.1f')
                    ytickformat('%.1f')
                    ax.FontName='Times';
                    ax.XAxis.FontSize = 16.5;
                    ax.YAxis.FontSize = 16.5;
                    ax.LineWidth = 1.5;


                    C = [pwd,"/",output,"/","Rog_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".png"];
                    str = join(C,"")
                    exportgraphics(ax ,str,'Resolution',125)
                    hold off
                    RR=Newt;
                    RR(:,2)=Rg2x(1:SubValue(L));
                    C = [pwd,"/",output,"/","Rog_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".csv"];
                    str = join(C,"")
                    writematrix(RR,str)



                    plot(Newt,(Rg2x(1:SubValue(L))),'color','black','LineWidth',2)
                    xlabel('time, s'), ylabel('<R_g^2>, Å^2')
                    set(gca,'FontWeight','bold')
                    ax = gca;
                    grid(ax,'on')
                    set(gca,'GridLineStyle','  -  -')
                    xtickformat('%.1f')
                    ytickformat('%.1f')
                    ax.FontName='Times';
                    ax.XAxis.FontSize = 16.5;
                    ax.YAxis.FontSize = 16.5;
                    ax.LineWidth = 1.5;



                    hold off
                    plot(Newt,(Tort2x(1:SubValue(L))),'color','black','LineWidth',2)
                    xlabel('time, s'), ylabel('<τ>, Å')
                    set(gca,'FontWeight','bold')
                    ax = gca;
                    grid(ax,'on')
                    set(gca,'GridLineStyle','  -  -')
                    xtickformat('%.1f')
                    ytickformat('%.1f')
                    ax.FontName='Times';
                    ax.XAxis.FontSize = 16.5;
                    ax.YAxis.FontSize = 16.5;
                    ax.LineWidth = 1.5;


                    C = [pwd,"/",output,"/","Tort_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),".png"];
                    str = join(C,"")
                    exportgraphics(ax ,str,'Resolution',125)
                    hold off

                    RRR=Newt;
                    RRR(:,2)=Tort2x(1:SubValue(L));
                    C = [pwd,"/",output,"/","Tort_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),".csv"];
                    str = join(C,"")
                    writematrix(RRR,str)

                    hold off
                    plot(Newt,(Rgperpx(1:SubValue(L))),'color','black','LineWidth',2)
                    xlabel('time, s'), ylabel('<R_g_\perp^2>, Å^2')
                    set(gca,'FontWeight','bold')
                    ax = gca;
                    grid(ax,'on')
                    set(gca,'GridLineStyle','  -  -')
                    xtickformat('%.1f')
                    ytickformat('%.1f')
                    ax.FontName='Times';
                    ax.XAxis.FontSize = 16.5;
                    ax.YAxis.FontSize = 16.5;
                    ax.LineWidth = 1.5;


                    C = [pwd,"/",output,"/","Rogperp_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),".png"];
                    str = join(C,"")
                    exportgraphics(ax ,str,'Resolution',125)
                    hold off

                    RRR=Newt;
                    RRR(:,2)=Rgperpx(1:SubValue(L));
                    C = [pwd,"/",output,"/","Rogperp_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),".csv"];
                    str = join(C,"")
                    writematrix(RRR,str)

                    plot(Newt(1:end-1),diff(Rgperpx(1:SubValue(L))),'color','black','LineWidth',2)
                    xlabel('time, s'), ylabel('d<R_g_\perp^2>/dt, Å^2')
                    set(gca,'FontWeight','bold')
                    ax = gca;
                    grid(ax,'on')
                    set(gca,'GridLineStyle','  -  -')
                    xtickformat('%.1f')
                    ytickformat('%.1f')
                    ax.FontName='Times';
                    ax.XAxis.FontSize = 16.5;
                    ax.YAxis.FontSize = 16.5;
                    ax.LineWidth = 1.5;


                    C = [pwd,"/",output,"/","dRog_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),".png"];
                    str = join(C,"")
                    exportgraphics(ax ,str,'Resolution',125)

                    RR=Newt;
                    RR(:,2)=Rgx(1:SubValue(L))';
                    C = [pwd,"/",output,"/","dRog_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),".csv"];
                    str = join(C,"")
                    writematrix(RR,str)

                    hold off
                    close(gcf)
                    plot(Newt,(Rgzx(1:SubValue(L))),'color','black','LineWidth',2)
                    xlabel('time, s'), ylabel('<R_g_z^2>, Å^2')
                    set(gca,'FontWeight','bold')
                    ax = gca;
                    grid(ax,'on')
                    set(gca,'GridLineStyle','  -  -')
                    xtickformat('%.1f')
                    ytickformat('%.1f')
                    ax.FontName='Times';
                    ax.XAxis.FontSize = 16.5;
                    ax.YAxis.FontSize = 16.5;
                    ax.LineWidth = 1.5;


                    C = [pwd,"/",output,"/","Rogz_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),".png"];
                    str = join(C,"")
                    exportgraphics(ax ,str,'Resolution',125)
                    hold off
                    close(gcf)
                    RRR=Newt;
                    RRR(:,2)=Rgzx(1:SubValue(L));
                    C = [pwd,"/",output,"/","Rogz_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),".csv"];
                    str = join(C,"")
                    writematrix(RRR,str)


                    SubSubValue =1;
                    x=0;
                    for xx=1:reducingindex:SubValue(L)
                        x=x+1;
                    end
                    SubSubValue =x;

                    close(gcf);
                    xlabel('Depth, Å'), ylabel('Polymer Brush Density, <ϕ_p>')
                    set(gca,'FontWeight','bold')
                    ax = gca;
                    grid(ax,'on')
                    set(gca,'GridLineStyle','  --')
                    xtickformat('%.2f')
                    ytickformat('%.2f')
                    ax.FontName='Times';
                    ax.XAxis.FontSize = 16.5;
                    ax.YAxis.FontSize = 16.5;
                    ax.LineWidth = 1.5;


                    N=length(DPx(:,1));
                    %xx=(1:length(DPx(1,:)))*2;
                    PPP=turbo(N);
                    %     DPmax=8;
                    hold on

                    maxstart=max(DPx(1,:));
                    maxend=max(DPx(SubSubValue,:));

                    for i=1:length(DPx(1:SubSubValue,1))
                        xx=(1:length(DPx(1,:)))*2.35;
                        plhdr = movmean(DPx(i,1:end),2);

                        y=plhdr(1:end);
                        x=xx(1:end);
                        start = 1;
                        stop =16;
                        width=5;
                        y(x>start & x<stop)=[];
                        x(x>start & x<stop)=[];
                        start2 = 102;
                        stop2 = x(end);
                        y(x>start2 & x<stop2)=[];
                        x(x>start2 & x<stop2)=[];


                        % map to new xaxis, leaving a space 'width' wide

                        x2=x;
                        x2(x2>=stop)=x2(x2>=stop)-(stop-start-width);
                        x2(end)=x2(end-1)+x2(end-2)-x2(end-3);
                        h=plot(x2,y,'LineWidth',2,'Color',PPP(i,:));
                        if((maxstart>=maxend && i==1 )||(maxend>maxstart && i==length(DPx(1:SubSubValue,1))))
                            ytick=get(gca,'YTick');


                            t11=text(.75+start+width/2,ytick(1),'▬','fontsize',8,'Color','white');
                            t22=text(.75+start+width/2,ytick(max(length(ytick))),'▬','fontsize',8,'Color','white');

                            t1=text(start+width/2,ytick(1),'//','fontsize',18);
                            t2=text(start+width/2,ytick(max(length(ytick))),'//','fontsize',18);

                            t11=text(.75+3+x2(end-1),ytick(1),'▬','fontsize',8,'Color','white');
                            t22=text(.75+3+x2(end-1),ytick(max(length(ytick))),'▬','fontsize',8,'Color','white');

                            t1=text(3+x2(end-1),ytick(1),'//','fontsize',18);
                            t2=text(3+x2(end-1),ytick(max(length(ytick))),'//','fontsize',18);
                            % For y-axis breaks, use set(t1,'rotation',270);
                            % remap tick marks, and 'erase' them in the gap
                            xtick=get(gca,'XTick');
                            dtick=xtick(2)-xtick(1);
                            gap=floor(width/dtick);
                            last=max(xtick(xtick<=start));          % last tick mark in LH dataset
                            next=min(xtick(xtick>=(last+dtick*(1+gap))));   % first tick mark within RH dataset
                            offset=size(x2(x2>last&x2<next),2)*(x(2)-x(1));
                            for i=1:sum(xtick>(last+gap))
                                xtick(find(xtick==last)+i+gap)=stop+offset+dtick*(i-1);
                            end
                            xtick(end)=x(end);
                            for i=1:length(xtick)
                                if xtick(i)>last&xtick(i)<next
                                    xticklabel{i}=sprintf('%d',[]);
                                else
                                    xticklabel{i}=sprintf('%d',xtick(i));
                                end
                            end
                            set(gca,'xticklabel',xticklabel);
                        end
                    end
                    %plot([0,8],[.06,.06],'b:','LineWidth',8)
                    hold off


                    C = [pwd,"/",output,"/","Density_Poly_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".png"];
                    str = join(C,"")
                    exportgraphics(ax ,str,'Resolution',125)

                    hold off



                    hold off
                    close(gcf)
                    xlabel('Depth, Å'), ylabel('Solvent Density , <ϕ_s>')
                    set(gca,'FontWeight','bold')
                    ax = gca;
                    grid(ax,'on')
                    set(gca,'GridLineStyle','  -  -')
                    xtickformat('%.2f')
                    ytickformat('%.2f')
                    ax.FontName='Times';
                    ax.XAxis.FontSize = 16.5;
                    ax.YAxis.FontSize = 16.5;
                    ax.LineWidth = 1.5;

                    hold on
                    %     plot(movmean(DPSx(1,:),2))
                    %     for i=1:length(DPx(:,1))
                    %         plot(movmean(DPSx(i,:),2))
                    %     end
                    %
                    N=length(DPSx(:,1));
                    %xx=(1:length(DPSx(1,:)))*2;
                    PPP=turbo(N);
                    %DPmax=4;
                    hold on

                    maxstart=max(DPSx(1,:));
                    maxend=max(DPSx(SubSubValue,:));

                    for i=1:length(DPSx(1:SubSubValue,1))
                        %         plhdrDP=movmean(DPSx(i,:),2);
                        %         if max(plhdrDP)~=0
                        %             plot([0,xx(1,DPmax)],[plhdrDP(1,DPmax),plhdrDP(1,DPmax)],':','LineWidth',2,'Color',PPP(i,:))
                        %             %plot(xx(1,1:8),plhdrDP(1,1:8))
                        %             plot(xx(1,DPmax:69),plhdrDP(1,DPmax:69),'LineWidth',2,'Color',PPP(i,:))
                        %             plot([xx(1,69),xx(1,end)],[plhdrDP(1,69),plhdrDP(1,69)],':','LineWidth',2,'Color',PPP(i,:))
                        %             i
                        %         end

                        xx=(1:length(DPSx(i,:)))*2.35;
                        plhdr = movmean(DPSx(i,1:end),2);

                        y=plhdr(1:end);
                        x=xx(1:end);
                        start = 0;
                        stop =8;
                        width=6;
                        %stop =10;
                        %width=1;
                        y(x>start & x<stop)=[];
                        x(x>start & x<stop)=[];
                        start2 = 102;
                        stop2 = x(end)+1;
                        y(x>start2 & x<stop2)=[];
                        x(x>start2 & x<stop2)=[];


                        % map to new xaxis, leaving a space 'width' wide
                        x2=x;
                        x2(x2>=stop)=x2(x2>=stop)-(stop-start-width);
                        x2(end)=x2(end-1)+x2(end-2)-x2(end-3);
                        h=plot(x2,y,'LineWidth',2,'Color',PPP(i,:));
                        if((maxstart>=maxend && i==1 )||(maxend>maxstart && i==length(DPSx(1:SubSubValue,1))))
                            ytick=get(gca,'YTick');

                            t11=text(1.2+start+width/2,ytick(1),'▬','fontsize',8,'Color','white');
                            t22=text(1.2+start+width/2,ytick(max(length(ytick))),'▬','fontsize',8,'Color','white');

                            t1=text(start+width/2,ytick(1),'//','fontsize',18);
                            t2=text(start+width/2,ytick(max(length(ytick))),'//','fontsize',18);

                            t11=text(1.2+x2(end),ytick(1),'▬','fontsize',8,'Color','white');
                            t22=text(1.2+x2(end),ytick(max(length(ytick))),'▬','fontsize',8,'Color','white');

                            t1=text(x2(end),ytick(1),'//','fontsize',18);
                            t2=text(x2(end),ytick(max(length(ytick))),'//','fontsize',18);
                            % For y-axis breaks, use set(t1,'rotation',270);
                            % remap tick marks, and 'erase' them in the gap
                            xtick=get(gca,'XTick');
                            dtick=xtick(2)-xtick(1);
                            gap=floor(width/dtick);
                            last=max(xtick(xtick<=start));          % last tick mark in LH dataset
                            next=min(xtick(xtick>=(last+dtick*(1+gap))));   % first tick mark within RH dataset
                            offset=size(x2(x2>last&x2<next),2)*(x(2)-x(1));
                            for i=1:sum(xtick>(last+gap))
                                xtick(find(xtick==last)+i+gap)=stop+offset+dtick*(i-1);
                            end
                            xtick(end)=xx(end);
                            for i=1:length(xtick)
                                if xtick(i)>last&xtick(i)<next
                                    xticklabel{i}=sprintf('%d',[]);
                                else
                                    xticklabel{i}=sprintf('%d',xtick(i));
                                end
                            end
                            set(gca,'xticklabel',xticklabel);

                        end
                    end
                    hold off

                    C = [pwd,"/",output,"/","Density_Solvent_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".png"];
                    str = join(C,"")
                    exportgraphics(ax ,str,'Resolution',125)



                    plot(shortt(1:SubSubValue),(Heightx(1:SubSubValue))*2.35,'color','black','LineWidth',2)
                    xlabel('time, s'), ylabel('Brush Height, <h>, Å')
                    set(gca,'FontWeight','bold')
                    ax = gca;
                    grid(ax,'on')
                    set(gca,'GridLineStyle','  -  -')
                    xtickformat('%.2f')
                    ytickformat('%.2f')
                    ax.FontName='Times';
                    ax.XAxis.FontSize = 16.5;
                    ax.YAxis.FontSize = 16.5;
                    ax.LineWidth = 1.5;

                    ax.FontName='Times';
                    ax.XAxis.FontSize = 16;
                    ax.YAxis.FontSize = 16;


                    hold off
                    plot(shortt(1:SubSubValue),(Heightx(1:SubSubValue))*2.35,'color','black','LineWidth',2)
                    xlabel('time, s'), ylabel('Expected Brush Height, <h>, Å')
                    set(gca,'FontWeight','bold')
                    ax = gca;
                    grid(ax,'on')
                    set(gca,'GridLineStyle','  -  -')
                    xtickformat('%.2f')
                    ytickformat('%.2f')
                    ax.FontName='Times';
                    ax.XAxis.FontSize = 16.5;
                    ax.YAxis.FontSize = 16.5;
                    ax.LineWidth = 1.5;


                    C = [pwd,"/",output,"/","Height_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".png"];
                    str = join(C,"")
                    exportgraphics(ax ,str,'Resolution',125)
                    hold off

                    RRR=shortt(1:SubSubValue);
                    RRR(:,2)=Heightx(1:SubSubValue);
                    C = [pwd,"/",output,"/","Height_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".csv"];
                    str = join(C,"")
                    writematrix(RRR,str)

                    hold off
                    plot(shortt(1:SubSubValue),(Widthx(1:SubSubValue)),'color','black','LineWidth',2)
                    xlabel('time, s'), ylabel('Expected Interfacial Width, <W>,Å*dϕ_p/dÅ')
                    set(gca,'FontWeight','bold')
                    ax = gca;
                    grid(ax,'on')
                    xtickformat('%.2f')
                    ytickformat('%.2f')
                    ax.FontName='Times';
                    ax.XAxis.FontSize = 16;
                    ax.YAxis.FontSize = 16;

                    C = [pwd,"/",output,"/","Width_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".png"];
                    str = join(C,"")
                    exportgraphics(ax ,str,'Resolution',125)
                    hold off

                    RRR=shortt(1:SubSubValue);
                    RRR(:,2)=Widthx(1:SubSubValue)';
                    C = [pwd,"/",output,"/","Width_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".csv"];
                    str = join(C,"")
                    writematrix(RRR,str)

                    hold off
                    plot(shortt(1:SubSubValue),(Abx(1:SubSubValue))*100,'color','black','LineWidth',2)
                    xlabel('time, s'), ylabel('Expected Absorbance, <%_a_b>')
                    set(gca,'FontWeight','bold')
                    ax = gca;
                    grid(ax,'on')
                    xtickformat('%.2f')
                    ytickformat('%.2f')
                    ax.FontName='Times';
                    ax.XAxis.FontSize = 16;
                    ax.YAxis.FontSize = 16;

                    C = [pwd,"/",output,"/","Absorbance_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".png"];
                    str = join(C,"")
                    exportgraphics(ax ,str,'Resolution',125)
                    hold off

                    RRR=shortt(1:SubSubValue);
                    RRR(:,2)=Abx(1:SubSubValue);
                    C = [pwd,"/",output,"/","Absorbance_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".csv"];
                    str = join(C,"")
                    writematrix(RRR,str)


                    hold off
                    plot(shortt(1:SubSubValue),(Adx(1:SubSubValue))*100,'color','black','LineWidth',2)
                    xlabel('time, s'), ylabel('Expected Adsorbance, <%_a_d>')
                    set(gca,'FontWeight','bold')
                    ax = gca;
                    grid(ax,'on')
                    xtickformat('%.2f')
                    ytickformat('%.1f')
                    ax.FontName='Times';
                    ax.XAxis.FontSize = 16;
                    ax.YAxis.FontSize = 16;

                    C = [pwd,"/",output,"/","Adsorbance_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".png"];
                    str = join(C,"")
                    exportgraphics(ax ,str,'Resolution',125)
                    hold off

                    RRR=shortt(1:SubSubValue);
                    RRR(:,2)=Adx(1:SubSubValue);
                    C = [pwd,"/",output,"/","Adsorbance_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".csv"];
                    str = join(C,"")
                    writematrix(RRR,str)

                    hold off
                    plot(shortt(1:SubSubValue),(Addx(1:SubSubValue)),'color','black','LineWidth',2)
                    xlabel('time, s'), ylabel('Expected Adsorbance, <%_a_d> Å^2')
                    set(gca,'FontWeight','bold')
                    ax = gca;
                    grid(ax,'on')
                    xtickformat('%.2f')
                    ytickformat('%.1f')
                    ax.FontName='Times';
                    ax.XAxis.FontSize = 16;
                    ax.YAxis.FontSize = 16;

                    C = [pwd,"/",output,"/","Also_Adsorbance_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".png"];
                    str = join(C,"")
                    exportgraphics(ax ,str,'Resolution',125)
                    hold off

                    RRR=shortt(1:SubSubValue);
                    RRR(:,2)=Addx(1:SubSubValue);
                    C = [pwd,"/",output,"/","Also_Adsorbance_CycloBrush_",string(realV),"_T_",string(Temp),"_sub_",string(L),"_",string(VV),".csv"];
                    str = join(C,"")
                    writematrix(RRR,str)
                    hold off

                end
            end

            close(gcf)

            RepMicroHelp = [0 0 0 0];
            for i=1:reducingindex:length(P(:,1))
                RepMicroHelp(i,:) = [Ex(i)-120000 Sx(i) Rgx(i) Tortx(i)];
            end
            C = [pwd,"/",output,"/","RepMicroHelp_",string(realV),"_T_",string(Temp),"_",string(VV),".txt"];
            str = join(C,"")
            writematrix(RepMicroHelp,str)
            hold off


            RepMicroHelp = zeros(length(DPx),length(DPx(1,:)));
            counter=1;
            for i=1:reducingindex:length(P(:,1))
                RepMicroHelp(i,:) = DPx(counter,:);
                counter=counter+1;
            end
            C = [pwd,"/",output,"/","RepDPMicroHelp_",string(realV),"_T_",string(Temp),"_",string(VV),".txt"];
            str = join(C,"")
            writematrix(RepMicroHelp,str)

            RepMicroHelp = zeros(length(DPx),length(DPx(1,:)));
            counter=1;
            for i=1:reducingindex:length(P(:,1))
                RepMicroHelp(i,:) = DPSx(counter,:);
                counter=counter+1;
            end
            C = [pwd,"/",output,"/","RepDPSMicroHelp_",string(realV),"_T_",string(Temp),"_",string(VV),".txt"];
            str = join(C,"")
            writematrix(RepMicroHelp,str)

            hold off

            if(Temp<=10)
                WW=zeros(length(P(:,1)),2);
                WW(1,1)=(Ex(1));
                WW(1,2)=(Rgx(1));
                WW(1,3)=(Tortx(1));

                WWSE=zeros(4,2);
                WWSE(1,1)=(Sx(1));
                WWSE(1,2)=(Ex(1));

                if(ceil(Ex(1))>=Ex(1))
                    rrr=ceil(Ex(1));
                end

                iii=2;
                for UU=1:length(Ex(1,:))
                    if(rrr>=Ex(UU))
                        WW(iii,1)=(round(Ex(UU)));
                        WW(iii,2)=(Rgx(UU));
                        WW(iii,3)=(Tortx(UU));

                        WWSE(iii,1)=(Sx(UU));
                        WWSE(iii,2)=(Ex(UU));

                        Ex(UU)

                        iii=iii+1;
                        rrr=round(Ex(UU))-1;
                    end
                end

                WW(iii,1)=(Ex(length(Ex(1,:))));
                WW(iii,2)=(Tortx(length(Ex(1,:))));
                WW(iii,3)=(Rgx(length(Ex(1,:))));
                WWSE(iii,1)=(Sx(length(Ex(1,:))));
                WWSE(iii,2)=(Ex(length(Ex(1,:))));

                plot(ExvS(:,1),ExvS(:,2),'black','LineWidth',2)
                hold on
                plot(Sx,Ex,'color','red','LineWidth',2)
                hold on
                scatter(WWSE(1,1),WWSE(1,2),5,'MarkerEdgeColor',[.64 .64 .64],'MarkerFaceColor',[.64 .64 .64],'LineWidth',1)
                scatter(WWSE(end,1),WWSE(end,2),5,'MarkerEdgeColor',[.64 .64 .64],'MarkerFaceColor',[.64 .64 .64],'LineWidth',1)
                ylabel('Energy, J/mol'), xlabel('Entropy, J per mol/K')
                set(gca,'FontWeight','bold')
                ax = gca;
                grid(ax,'on')
                xtickformat('%.1f')
                ytickformat('%.1f')
                ax.FontName='Times';
                ax.XAxis.FontSize = 16;
                ax.YAxis.FontSize = 16;
                hold off
                C = [pwd,"/",output,"/","ExvSx_Marked_CycloBrush_",string(realV),"_T_",string(Temp),".png"];
                str = join(C,"")
                exportgraphics(ax ,str,'Resolution',125)
                Temp
                %pause(10)
                C = [pwd,"/",output,"/","Combined_CycloBrush_",string(realV),"_T_",string(Temp),".txt"];
                str = join(C,"")
                writematrix(WW,str)
                RRR=transpose(Sx);
                RRR(:,2)=Ex;
                C = [pwd,"/",output,"/","ExvSx_S_CycloBrush_",string(realV),"_T_",string(Temp),".txt"];
                str = join(C,"")
                writematrix(RRR,str)
            end

            if(Temp>=100)
                WW=zeros(length(P(:,1)),2);
                WW(1,1)=(Ex(1));
                WW(1,2)=(Rgx(1));
                WW(1,3)=(Tortx(1));

                WWSE=zeros(4,2);
                WWSE(1,1)=(Sx(1));
                WWSE(1,2)=(Ex(1));

                if(ceil(Ex(1))>=Ex(1))
                    rrr=ceil(Ex(1));
                end
                iii=2;
                for UU=1:length(Ex(1,:))
                    if(rrr<=Ex(UU))
                        WW(iii,1)=(round(Ex(UU)));
                        WW(iii,2)=(Rgx(UU));
                        WW(iii,3)=(Tortx(UU));

                        WWSE(iii,1)=(Sx(UU));
                        WWSE(iii,2)=(Ex(UU));

                        Ex(UU)

                        iii=iii+1;
                        rrr=round(Ex(UU))+1;
                    end
                end
                if(iii==2)
                    for UU=1:length(Ex(1,:))
                        if(rrr>=Ex(UU))
                            WW(iii,1)=(round(Ex(UU)));
                            WW(iii,2)=(Rgx(UU));
                            WW(iii,3)=(Tortx(UU));

                            WWSE(iii,1)=(Sx(UU));
                            WWSE(iii,2)=(Ex(UU));

                            Ex(UU)

                            iii=iii+1;
                            rrr=round(Ex(UU))+1;
                        end
                    end
                end
                WW(iii,1)=(Ex(length(Ex(1,:))));
                WW(iii,2)=(Tortx(length(Ex(1,:))));
                WW(iii,3)=(Rgx(length(Ex(1,:))));
                WWSE(iii,1)=(Sx(length(Ex(1,:))));
                WWSE(iii,2)=(Ex(length(Ex(1,:))));

                plot(ExvS(:,1),ExvS(:,2),'black','LineWidth',2)
                hold on
                plot(Sx,Ex,'color','red','LineWidth',2)
                hold on
                scatter(WWSE(1,1),WWSE(1,2),15,'MarkerEdgeColor',[.64 .64 .64],'MarkerFaceColor',[.64 .64 .64],'LineWidth',1)
                scatter(WWSE(end,1),WWSE(end,2),15,'MarkerEdgeColor',[.64 .64 .64],'MarkerFaceColor',[.64 .64 .64],'LineWidth',1)
                ylabel('Energy, J/mol'), xlabel('Entropy, J per mol/K')
                set(gca,'FontWeight','bold')
                ax = gca;
                grid(ax,'on')
                xtickformat('%.1f')
                ytickformat('%.1f')
                ax.FontName='Times';
                ax.XAxis.FontSize = 16;
                ax.YAxis.FontSize = 16;
                hold off
                C = [pwd,"/",output,"/","ExvSx_Marked_CycloBrush_",string(realV),"_T_",string(Temp),".png"];
                str = join(C,"")
                exportgraphics(ax ,str,'Resolution',125)
                Temp
                %pause(10)
                C = [pwd,"/",output,"/","Combined_CycloBrush_",string(realV),"_T_",string(Temp),".txt"];
                str = join(C,"")
                writematrix(WW,str)
                RRR=transpose(Sx);
                RRR(:,2)=Ex;
                C = [pwd,"/",output,"/","ExvSx_S_CycloBrush_",string(realV),"_T_",string(Temp),".txt"];
                str = join(C,"")
                writematrix(RRR,str)
            end
            %RR=t;
            %RR(:,2)=W';
            %writematrix(RR,"Tort_CycloBrush_Data_T_25.csv")
        end
        %end
    end
    %{
Beta=Sx./Ex;

Aa=zeros(1,length(P(:,1)));
Es=zeros(1,length(P(:,1)));
E2s=zeros(1,length(P(:,1)));
for U=1:length(P(:,1))
    for UU=1:length(P(1,:))
        Es(U)=Es(U)+P(U,UU)*E(UU);
        E2s(U)=E2s(U)+P(U,UU)*(E(UU))^2;
    end
     Aa(U)=E2s(U)-Es(U)^2;
end


Res_Aa=0;
Res_Es=0;
Res_E2s=0;
for UU=1:length(Pse(1,:))
    Res_Es=Res_Es+Pse(1,UU)*E(UU);
    Res_E2s=Res_E2s+Pse(1,UU)*(E(UU))^2;
end
Res_Aa=Res_E2s-Res_Es^2;
    %}
    %% Finding
    %{
J=importdata("34x34_R_Levels.xlsx");
Surface = zeros(6,length(Pse));
Combined = zeros(12,36);
Combined_Further = zeros(2,36);
CombinedIndex = 0;
difference=10000;
increment = 1;
incrementc=1;

for c=1:length(Pse)
        Surface(1,c)=10000000;
        Surface(2,c)=10000000;
        Surface(3,c)=10000000;
        Surface(4,c)=10000000;
        Surface(5,c)=10000000;
        Surface(6,c)=10000000;
end

for b=1:1:length(t)
    
    for c=1:length(Pse)
        Surface(1,c)=10000000;
        Surface(2,c)=10000000;
        Surface(3,c)=10000000;
        Surface(4,c)=10000000;
        Surface(5,c)=10000000;
        Surface(6,c)=10000000;
    end
    
    for c=1:length(Pse)
        if((WW(b)>(SK(c)-SK(c)*.05)) && (WW(b)<(SK(c)+SK(c)*.05)))
            %if((WWW(b)>(SKU(c)-SKU(c)*.25)) && (WWW(b)<(SKU(c)+SKU(c)*.25)))
                if((WT(b)>(T(c)-T(c)*.1)) && (WT(b)<(T(c)+T(c)*.1 )))
            Surface(1,increment)=SK(c);
            Surface(2,increment)=SKU(c);
            Surface(3,increment)=P(c);
            Surface(4,increment)=c;
            Surface(5,increment)=SQT(c);
            Surface(6,increment)=T(c);
            increment = increment +1;
                end
            %end
        end
    end
    store=1
    for a=1:length(Pse)
        if( difference>=(abs(Surface(2,a)-WWW(b))))
            difference = (abs(Surface(2,a)-WWW(b)));
            store=a;
        end
    end
    
    if(Surface(4,store)~=10000000)
    Combined(1,incrementc) = Surface(1,store);
    Combined(2,incrementc) = Surface(2,store);
    Combined(3,incrementc) = Surface(3,store);
    Combined(4,incrementc) = Surface(4,store);
    Combined(5,incrementc) = J(Surface(4,store));
    Combined(6,incrementc) = t(b);
    Combined(7,incrementc) = SQT(Surface(4,store));
    Combined(8,incrementc) = T(Surface(4,store));
    Combined(9,incrementc) = WW(b);
    Combined(10,incrementc) = WWW(b);
    Combined(11,incrementc) = W(b);
    Combined(12,incrementc) = WT(b);
    
    Combined_Further(1,incrementc)= J(Surface(4,store));
    Combined_Further(2,incrementc)= T(Surface(4,store));
    end

    difference=10000;
    increment =1;
    disp(b);
    incrementc=incrementc+1;
    
end

PPHH=Combined_Further(2,:).';
for q=1:length(Combined(11,:))

    CombinedIndex(q,1)=Combined(5,q);
    
end
    %}
    %%


    plot(ExvS(:,1)*scale*kB,ExvS(:,2)*scale,'Color',[.5 .5 .5],'LineWidth',1.5)
    hold on

    path=zeros(30000,2);
    d=winter(length(c_info(:,1))+1);
    for I=2:reduce:(length(c_info(:,1))-1)
        %for I=1:1
        dummy = importdata(strcat(pwd,"/",string(VV),"/",string(I),"/","ExvSx_CycloBrush_",num2str(I),"_T_300.csv"));
        for II = 1:length(path(:,1))
            if II<=length(dummy(:,1))
                path(II,2*(I+1)-1) = dummy(II,1);
                path(II,2*(I+1)) = dummy(II,2);
            else
                path(II,2*(I+1)-1) = dummy(length(dummy(:,1)),1);
                path(II,2*(I+1)) = dummy(length(dummy(:,1)),2);
            end
        end
        p=plot(path(:,2*(I+1)-1)*scale*kB,path(:,2*(I+1))*scale*kB,'Color',d((I+1),:),'LineWidth',2);
        p.LineStyle = '--';
        hold on
    end
    %path=zeros(3000,2);
    loc=1;
    for K=1:length(ExvS(:,1))
        if(round(100/VV,1)<=ExvS(K,3))
            loc =K;
            break;
        end
    end

    scatter(ExvS(loc,1)*scale*kB,ExvS(loc,2)*scale,'red','LineWidth',8)

    % path1 = importdata("ExvSx_S_CycloBrush_"<>num2str(I)<>"_T_100.csv");
    % path2 = importdata("ExvSx_S_CycloBrush_2_T_100.csv");
    % path3 = importdata("ExvSx_S_CycloBrush_3_T_100.csv");
    % path4 = importdata("ExvSx_S_CycloBrush_4_T_100.csv");
    % path5 = importdata("ExvSx_S_CycloBrush_5_T_100.csv");
    % path6 = importdata("ExvSx_S_CycloBrush_6_T_100.csv");
    % path7 = importdata("ExvSx_S_CycloBrush_7_T_100.csv");
    % path8 = importdata("ExvSx_S_CycloBrush_8_T_100.csv");
    % path9 = importdata("ExvSx_S_CycloBrush_9_T_100.csv");
    % path10 = importdata("ExvSx_S_CycloBrush_10_T_100.csv");
    % path11 = importdata("ExvSx_S_CycloBrush_11_T_100.csv");
    % path12 = importdata("ExvSx_S_CycloBrush_12_T_100.csv");
    % path13 = importdata("ExvSx_S_CycloBrush_13_T_100.csv");
    % path14 = importdata("ExvSx_S_CycloBrush_14_T_100.csv");
    % path15 = importdata("ExvSx_S_CycloBrush_15_T_100.csv");
    % path16 = importdata("ExvSx_S_CycloBrush_16_T_100.csv");
    % path17 = importdata("ExvSx_S_CycloBrush_17_T_100.csv");
    % path18 = importdata("ExvSx_S_CycloBrush_18_T_100.csv");

    % repm1 = importdata("Combined_CycloBrush_1_T_100.txt");
    % repm2 = importdata("Combined_CycloBrush_9_T_100.txt");
    % repm3 = importdata("Combined_CycloBrush_17_T_5.txt");
    %
    % repmout1=[0 0 0 ; 0 0 0];
    % repmout2=[0 0 0 ; 0 0 0];
    % repmout3=[0 0 0 ;0 0 0];

    %     plot(ExvS(:,1)*scale,ExvS(:,2)*scale,'black','LineWidth',4)
    %     hold on
    % plot(path1(:,1)*scale,path1(:,2)*scale,'red','LineWidth',4)
    % hold on
    % plot(path2(:,1)*scale,path2(:,2)*scale,'blue','LineWidth',2)
    % hold on
    % plot(path3(:,1)*scale,path3(:,2)*scale,'magenta','LineWidth',4)
    % hold on
    % plot(path4(:,1)*scale,path4(:,2)*scale,'cyan','LineWidth',2)
    % hold on
    % plot(path5(:,1)*scale,path5(:,2)*scale,'blue','LineWidth',2)
    % hold on
    % plot(path6(:,1)*scale,path6(:,2)*scale,'black','LineWidth',2)
    % hold on
    % plot(path7(:,1)*scale,path7(:,2)*scale,'red','LineWidth',4)
    % hold on
    % plot(path8(:,1)*scale,path8(:,2)*scale,'blue','LineWidth',2)
    % hold on
    % plot(path9(:,1)*scale,path9(:,2)*scale,'cyan','LineWidth',4)
    % hold on
    % plot(path10(:,1)*scale,path10(:,2)*scale,'Color',[255/255 165/255 0],'LineWidth',2)
    % hold on
    % plot(path11(:,1)*scale,path11(:,2)*scale,'cyan','LineWidth',2)
    % hold on
    % plot(path12(:,1)*scale,path12(:,2)*scale,'magenta','LineWidth',2)
    % hold on
    % plot(path13(:,1)*scale,path13(:,2)*scale,'cyan','LineWidth',2)
    % hold on
    % plot(path14(:,1)*scale,path14(:,2)*scale,'magenta','LineWidth',2)
    % hold on
    % plot(path15(:,1)*scale,path15(:,2)*scale,'cyan','LineWidth',2)
    % hold on
    % plot(path16(:,1)*scale,path16(:,2)*scale,'magenta','LineWidth',2)
    % hold on
    % plot(path17(:,1)*scale,path17(:,2)*scale,'Color',[255/255 165/255 0],'LineWidth',4)
    % hold on
    % plot(path18(:,1)*scale,path18(:,2)*scale,'cyan','LineWidth',2)
    % hold on
    % count=0;
    % for i=1:1:length(repm1(:,1))
    %     ccc=[3595859.85197660,75207105.9179414;
    %         3634297.45826208,82913505.0100918;
    %         3671844.36003283,90490688.7310587;
    %         3709394.12905239,98087410.8241246;
    %         3735896.36592203,103460368.408371;
    %         3763230.19215532,109018498.572579];
    %     if(count == 0 && repm1(i,1)>=ccc(count+1,2)/scale)
    %         repmout1(count+1,:) = repm1(i,:);
    %         count = count+1
    %     end
    %     if(count == 1&& repm1(i,1)>=ccc(count+1,2)/scale)
    %         repmout1(count+1,:) = repm1(i,:)
    %         count = count+1
    %     end
    %     if(count == 2 && repm1(i,1)>=ccc(count+1,2)/scale)
    %         repmout1(count+1,:) = repm1(i,:)
    %         count = count+1
    %     end
    %     if(count == 3&& repm1(i,1)>=ccc(count+1,2)/scale)
    %         repmout1(count+1,:) = repm1(i,:)
    %         count = count+1
    %     end
    %     if(count == 4 && repm1(i,1)>=ccc(count+1,2)/scale)
    %         repmout1(count+1,:) = repm1(i,:)
    %         count = count+1
    %     end
    %     if(count == 5 && repm1(i,1)>=ccc(count+1,2)/scale)
    %         repmout1(count+1,:) = repm1(i,:)
    %         count = count+1
    %     end
    % end
    %
    % count=0;
    % for i=1:1:length(path1(:,1))
    %     ccc=[3595859.85197660,75207105.9179414;
    %         3634297.45826208,82913505.0100918;
    %         3671844.36003283,90490688.7310587;
    %         3709394.12905239,98087410.8241246;
    %         3735896.36592203,103460368.408371;
    %         3763230.19215532,109018498.572579];
    %     if(count == 0 && path1(i,1)*scale>=ccc(count+1,1) && path1(i,2)*scale>=ccc(count+1,2))
    %         scatter(path1(i,1)*scale,path1(i,2)*scale,20,'MarkerEdgeColor',[.44 .44 .44],'MarkerFaceColor',[.44 .44 .44],'LineWidth',1)
    %         count = count+1
    %     end
    %     if(count == 1 && path1(i,1)*scale>=ccc(count+1,1) && path1(i,2)*scale>=ccc(count+1,2))
    %         scatter(path1(i,1)*scale,path1(i,2)*scale,20,'MarkerEdgeColor',[.44 .44 .44],'MarkerFaceColor',[.44 .44 .44],'LineWidth',1)
    %         count = count+1
    %     end
    %     if(count == 2 && path1(i,1)*scale>=ccc(count+1,1) && path1(i,2)*scale>=ccc(count+1,2))
    %         scatter(path1(i,1)*scale,path1(i,2)*scale,20,'MarkerEdgeColor',[.44 .44 .44],'MarkerFaceColor',[.44 .44 .44],'LineWidth',1)
    %         count = count+1
    %     end
    %     if(count == 3 && path1(i,1)*scale>=ccc(count+1,1) && path1(i,2)*scale>=ccc(count+1,2))
    %         scatter(path1(i,1)*scale,path1(i,2)*scale,20,'MarkerEdgeColor',[.44 .44 .44],'MarkerFaceColor',[.44 .44 .44],'LineWidth',1)
    %         count = count+1
    %     end
    %     if(count == 4 && path1(i,1)*scale>=ccc(count+1,1) && path1(i,2)*scale>=ccc(count+1,2))
    %         scatter(path1(i,1)*scale,path1(i,2)*scale,20,'MarkerEdgeColor',[.44 .44 .44],'MarkerFaceColor',[.44 .44 .44],'LineWidth',1)
    %         count = count+1
    %     end
    %     if(count == 5 && path1(i,1)*scale>=ccc(count+1,1) && path1(i,2)*scale>=ccc(count+1,2))
    %         scatter(path1(i,1)*scale,path1(i,2)*scale,20,'MarkerEdgeColor',[.44 .44 .44],'MarkerFaceColor',[.44 .44 .44],'LineWidth',1)
    %         count = count+1
    %     end
    % end
    %
    % hold on
    % count = 0;
    % for i=1:1:length(path9(:,1))
    %     if(count == 0 && path9(i,1)*scale>=2396836.48061096 && path9(i,2)*scale>=71804014.2144774)
    %         scatter(path9(i,1)*scale,path9(i,2)*scale,20,'MarkerEdgeColor',[.44 .44 .44],'MarkerFaceColor',[.44 .44 .44],'LineWidth',1)
    %         count = count+1
    %     end
    %     if(count == 1 && path9(i,1)*scale>=2877863.22786578 && path9(i,2)*scale>=103634854.418065)
    %         scatter(path9(i,1)*scale,path9(i,2)*scale,20,'MarkerEdgeColor',[.44 .44 .44],'MarkerFaceColor',[.44 .44 .44],'LineWidth',1)
    %         count = count+1
    %     end
    %     if(count == 2 && path9(i,1)*scale>=3332302.52456600 && path9(i,2)*scale>=133694043.752681)
    %         scatter(path9(i,1)*scale,path9(i,2)*scale,20,'MarkerEdgeColor',[.44 .44 .44],'MarkerFaceColor',[.44 .44 .44],'LineWidth',1)
    %         count = count+1
    %     end
    %     if(count == 3 && path9(i,1)*scale>=3755594.07231135 && path9(i,2)*scale<=161330327.985888)
    %         scatter(path9(i,1)*scale,path9(i,2)*scale,20,'MarkerEdgeColor',[.44 .44 .44],'MarkerFaceColor',[.44 .44 .44],'LineWidth',1)
    %         count = count+1
    %     end
    %     if(count == 4 && path9(i,1)*scale>=3761842.27397531 && path9(i,2)*scale<=133388041.417777)
    %         scatter(path9(i,1)*scale,path9(i,2)*scale,20,'MarkerEdgeColor',[.44 .44 .44],'MarkerFaceColor',[.44 .44 .44],'LineWidth',1)
    %         count = count+1
    %     end
    %     if(count == 5 && path9(i,1)*scale>=3766430.14350796 && path9(i,2)*scale<=111169114.067502)
    %         scatter(path9(i,1)*scale,path9(i,2)*scale,20,'MarkerEdgeColor',[.44 .44 .44],'MarkerFaceColor',[.44 .44 .44],'LineWidth',1)
    %         count = count+1
    %     end
    %
    % end
    % hold on
    % count=0;
    % for i=1:1:length(path17(:,1))
    %     ccc=[3512000,73940000;
    %     2681423.17834142,68366373.0202542;
    %     1959564.19564376,70080598.3465202;
    %     2014815.08090325,88534915.2797944;
    %     2352572.11249788,110332921.298343;
    %     2755537.10542399,134246815.880808];
    %     ccc=flip(ccc);
    %
    %     if(count == 0 && path17(i,1)*scale<=ccc(count+1,1) && path17(i,2)*scale<=ccc(count+1,2))
    %         scatter(path17(i,1)*scale,path17(i,2)*scale,20,'MarkerEdgeColor',[.44 .44 .44],'MarkerFaceColor',[.44 .44 .44],'LineWidth',2)
    %         count = count+1
    %     end
    %     if(count == 1 && path17(i,1)*scale<=ccc(count+1,1) && path17(i,2)*scale<=ccc(count+1,2))
    %         scatter(path17(i,1)*scale,path17(i,2)*scale,20,'MarkerEdgeColor',[.44 .44 .44],'MarkerFaceColor',[.44 .44 .44],'LineWidth',2)
    %         count = count+1
    %     end
    %     if(count == 2 && path17(i,1)*scale<=ccc(count+1,1) && path17(i,2)*scale<=ccc(count+1,2))
    %         scatter(path17(i,1)*scale,path17(i,2)*scale,20,'MarkerEdgeColor',[.44 .44 .44],'MarkerFaceColor',[.44 .44 .44],'LineWidth',2)
    %         count = count+1
    %     end
    %     if(count == 3 && path17(i,1)*scale>=ccc(count+1,1) && path17(i,2)*scale<=ccc(count+1,2))
    %         scatter(path17(i,1)*scale,path17(i,2)*scale,20,'MarkerEdgeColor',[.44 .44 .44],'MarkerFaceColor',[.44 .44 .44],'LineWidth',2)
    %         count = count+1
    %     end
    %     if(count == 4 && path17(i,1)*scale>=ccc(count+1,1) && path17(i,2)*scale>=ccc(count+1,2))
    %         scatter(path17(i,1)*scale,path17(i,2)*scale,20,'MarkerEdgeColor',[.44 .44 .44],'MarkerFaceColor',[.44 .44 .44],'LineWidth',2)
    %         count = count+1
    %     end
    %     if(count == 5 && path17(i,1)*scale>=ccc(count+1,1) && path17(i,2)*scale>=ccc(count+1,2))
    %         scatter(path17(i,1)*scale,path17(i,2)*scale,20,'MarkerEdgeColor',[.44 .44 .44],'MarkerFaceColor',[.44 .44 .44],'LineWidth',2)
    %         count = count+1
    %     end
    % end
    % hold on
    % for i=1:100:length(path4(:,1))
    %     scatter(path4(i,1)*scale,path4(i,2)*scale,15,'MarkerEdgeColor',[.44 .44 .44],'MarkerFaceColor',[.44 .44 .44],'LineWidth',1)
    % end
    xlabel('<S>, J per mol/K'), ylabel('<E>, J/mol')
    set(gca,'FontWeight','bold')
    ax = gca;
    grid(ax,'on')
    set(gca,'GridLineStyle','  -  -')
    xtickformat('%.1f')
    ytickformat('%.1f')
    ax.FontName='Times';
    ax.XAxis.FontSize = 16.5;
    ax.YAxis.FontSize = 16.5;
    ax.LineWidth = 1.5;

    C = ["ExvSx_1_",string(VV),"Combined_CycloBrush_CycloBrush.png"];
    str = join(C,"")
    exportgraphics(ax ,str,'Resolution',125)

    close(gcf);


    plot(ExvS(:,1)*scale*kB+.1*10^6*kB,ExvS(:,2)*scale-.1*10^7,'black','LineWidth',4)
    hold on
    plot(ExvS(:,1)*scale*kB,ExvS(:,2)*scale,'Color',[.5 .5 .5],'LineWidth',.5)
    hold on
    path=zeros(30000,2);
    d=winter(length(c_info(:,1))+1);
    for I=0:reduce:(length(c_info(:,1))-1)
        %for I=1:1
        dummy = importdata(strcat(pwd,"/",string(VV),"/",string(I),"/","ExvSx_CycloBrush_",num2str(I),"_T_300.csv"));
        for II = 1:length(path(:,1))
            if II<=length(dummy(:,1))
                path(II,2*(I+1)-1) = dummy(II,1);
                path(II,2*(I+1)) = dummy(II,2);
            else
                path(II,2*(I+1)-1) = dummy(length(dummy(:,1)),1);
                path(II,2*(I+1)) = dummy(length(dummy(:,1)),2);
            end
        end
        p=plot(path(:,2*(I+1)-1)*scale*kB,path(:,2*(I+1))*scale,'Color',d((I+1),:),'LineWidth',2);
        p.LineStyle = '--';
        hold on
    end
    %path=zeros(3000,2);
    loc=1;
    for K=1:length(ExvS(:,1))
        if(round(100/VV,1)<=ExvS(K,3))
            loc =K;
            break;
        end
    end

    scatter(ExvS(loc,1)*scale*kB,ExvS(loc,2)*scale,'red','LineWidth',5)

    xlabel('<S>, J per mol/K'), ylabel('<E>, J/mol')
    set(gca,'FontWeight','bold')
    ax = gca;
    grid(ax,'on')
    set(gca,'GridLineStyle','  -  -')
    xtickformat('%.1f')
    ytickformat('%.1f')
    ax.FontName='Times';
    ax.XAxis.FontSize = 16.5;
    ax.YAxis.FontSize = 16.5;
    ax.LineWidth = 1.5;

    C = ["ExvSx_ext_1_",string(VV),"Combined_CycloBrush_CycloBrush.png"];
    str = join(C,"")
    exportgraphics(ax ,str,'Resolution',125)



    if(VV==1)
        %for P=1:2:19
        % VV=P;
        close(gcf);

        plot(ExvS(:,1)*scale*kB,ExvS(:,2)*scale,'Color',[.5 .5 .5],'LineWidth',1.5)
        hold on

        d=autumn(3*length(c_info(:,1)));
        for I=start3:reduce:(length(c_info(:,1))+length(c_info2(:,1))-1)
            dummy = importdata(strcat(pwd,"/",string(VV),"/",string(I),"/","ExvSx_CycloBrush_",num2str(I),"_T_10.csv"));
            for II = 1:length(path(:,1))
                if II<=length(dummy(:,1))
                    path(II,2*I-1) = dummy(II,1);
                    path(II,2*I) = dummy(II,2);
                else
                    path(II,2*I-1) = dummy(length(dummy(:,1)),1);
                    path(II,2*I) = dummy(length(dummy(:,1)),2);
                end
            end
            plot(path(:,2*I-1)*scale*kB,path(:,2*I)*scale,'red','LineWidth',3)
            %ax.LineStyle = '--';
            hold on
        end

        loc=1;
        for K=1:length(ExvS(:,1))
            if(round(5/VV,1)<=ExvS(K,3))
                loc =K;
                break;
            end
        end
        hold on
        scatter(ExvS(loc,1)*scale*kB,ExvS(loc,2)*scale,'red','LineWidth',6)
        hold on
        xlabel('<S>, J per mol/K'), ylabel('<E>, J/mol')
        set(gca,'FontWeight','bold')
        ax = gca;
        grid(ax,'on')
        set(gca,'GridLineStyle','  -  -')
        xtickformat('%.1f')
        ytickformat('%.1f')
        ax.FontName='Times';
        ax.XAxis.FontSize = 16.5;
        ax.YAxis.FontSize = 16.5;
        ax.LineWidth = 1.5;

        C = ["ExvSx_2_",string(VV),"Combined_seq.png"];
        str = join(C,"")
        exportgraphics(ax ,str,'Resolution',125)

        close(gcf)


        plot(ExvS(:,1)*scale*kB,ExvS(:,2)*scale,'black','LineWidth',.5)
        hold on
        plot(ExvS(:,1)*scale*kB+.1*10^6*kB,ExvS(:,2)*scale-.1*10^7,'Color',[.5 .5 ,.5],'LineWidth',4)
        hold on

        d=autumn(3*length(c_info(:,1)));
        for I=start3:reduce:(length(c_info(:,1))+length(c_info2(:,1))-1)
            dummy = importdata(strcat(pwd,"/",string(VV),"/",string(I),"/","ExvSx_CycloBrush_",num2str(I),"_T_10.csv"));
            for II = 1:length(path(:,1))
                if II<=length(dummy(:,1))
                    path(II,2*I-1) = dummy(II,1);
                    path(II,2*I) = dummy(II,2);
                else
                    path(II,2*I-1) = dummy(length(dummy(:,1)),1);
                    path(II,2*I) = dummy(length(dummy(:,1)),2);
                end
            end
            p=plot(path(:,2*I-1)*scale*kB,path(:,2*I)*scale,'Color',d(I,:),'LineWidth',3);
            p.LineStyle = '--';
            hold on
        end

        loc=1;
        for K=1:length(ExvS(:,1))
            if(round(5/VV,1)<=ExvS(K,3))
                loc =K;
                break;
            end
        end
        hold on
        scatter(ExvS(loc,1)*scale*kB,ExvS(loc,2)*scale,'red','LineWidth',6)
        hold on
        xlabel('<S>, J per mol/K'), ylabel('<E>, J/mol')
        set(gca,'FontWeight','bold')
        ax = gca;
        grid(ax,'on')
        set(gca,'GridLineStyle','  -  -')
        xtickformat('%.1f')
        ytickformat('%.1f')
        ax.FontName='Times';
        ax.XAxis.FontSize = 16.5;
        ax.YAxis.FontSize = 16.5;
        ax.LineWidth = 1.5;

        C = ["ExvSx_ext_2_",string(VV),"Combined_seq.png"];
        str = join(C,"")
        exportgraphics(ax ,str,'Resolution',125)
    end
    %end
end