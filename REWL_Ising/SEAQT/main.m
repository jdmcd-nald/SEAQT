clear all
format longG
close(gcf)
global number_of_state energy_value number_of_degenerate_energy tau Concentration kB
global number_of_degenerate_energy_ln energy_value_ln Concentration_ln inc increment OldP Oldt PP PPP  ln_number_of_degenerate_energy beta_reservoir
global vq1 iiii rt tempSe rs


kB = 1;
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

plot(E,LNg,'color','black','LineWidth',2')
hold on
scatter(E,LNg,150,'MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor','black','LineWidth',1.5)
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


plot(E,LNg,'color','black','LineWidth',2')
hold on
for I=1:length(g(:,1))
    scatter(E(I),LNg(I),150,'MarkerFaceColor',colorrewl(I,:),'MarkerEdgeColor','black','LineWidth',2)
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


output = "Ising";
mkdir(output);

realV=1;


clear P
close(gcf)
close all

RR=ExvS(:,1)*kB;
RR(:,2)=ExvS(:,2);
C = [pwd,"/",output,"/","ExvSx",".xlsx"];
str = join(C,"")
writematrix(RR,str)




Temp=5;
ft=5;

minI=7;
maxI=20;
minII=7;
maxII=20;
reducep =30;


tspan=[0:.01:2]


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

plot(ExvS(:,1)*kB,ExvS(:,2),'Color',[.5 .5 .5],'LineWidth',.5)
hold on
plot(Sx*kB,Ex,'color','red','LineWidth',3)
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
plot(ExvS(:,1)*kB,ExvS(:,2),'Color',[.5 .5 .5],'LineWidth',.5)
hold on
plot(Sx*kB,Ex,'color','red','LineWidth',3)
xlabel('Entropy, J per mol/K'), ylabel('Energy, J/mol')
set(gca,'FontWeight','bold')
ax = gca;
grid(ax,'on')
set(gca,'GridLineStyle','  -  -')
xtickformat('%.0f')
ytickformat('%.0f')
ax.FontName='Times';

%     axis([ 5  28  ...
%     -50  0])
close(gcf)

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

plot(t,Ex,'color','black','LineWidth',2)
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

plot(t(1:end-1),diff(Ex),'color','black','LineWidth',2)
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
RR(:,2)=Ex';
C = [pwd,"/",output,"/","Ex",string(realV),"_T_",string(Temp),".csv"];
str = join(C,"")
writematrix(RR,str)

RR=t(2:end);
RR(:,2)=diff(Ex)';
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
HC=HC;

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

plot(t,Sx*kB,'color','black','LineWidth',2)
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

plot(t,Sx*kB-Sx(1)*kB,'color','black','LineWidth',2)
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

plot(t(1:end-1),diff(Sx*kB-Sx(1)*kB),'color','black','LineWidth',2)
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
RR(:,2)=(Sx*kB-Sx(1)*kB)';
C = [pwd,"/",output,"/","SxA_CycloBrush_",string(realV),"_T_",string(Temp),".csv"];
str = join(C,"")
writematrix(RR,str)

RR=t(1:end-1);
RR(:,2)=(diff(Sx*kB-Sx(1)*kB))';
C = [pwd,"/",output,"/","dSx_CycloBrush_",string(realV),"_T_",string(Temp),".csv"];
str = join(C,"")
writematrix(RR,str)

plot(t(1:end),(Sx*kB-Sx(1)*kB)-(Ex-Ex(1))/Temp,'color','black','LineWidth',2)
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
RR(:,2)=((Sx*kB-Sx(1)*kB)-(Ex-Ex(1))/Temp)';
C = [pwd,"/",output,"/","SxP_CycloBrush_",string(realV),"_T_",string(Temp),".csv"];
str = join(C,"")
writematrix(RR,str)

plot(t(1:end-1),diff((Sx*kB-Sx(1)*kB)-(Ex-Ex(1))/Temp),'color','black','LineWidth',2)
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
RR(:,2)=diff((Sx*kB-Sx(1)*kB)-(Ex-Ex(1))/Temp)';
C = [pwd,"/",output,"/","SigmaP_CycloBrush_",string(realV),"_T_",string(Temp),".csv"];
str = join(C,"")
writematrix(RR,str)

plot(ExvS(:,1)*kB,ExvS(:,2),'black','LineWidth',2)
hold on
plot(Sx*kB,Ex,'color','red','LineWidth',2)
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


hold off
C = [pwd,"/",output,"/","ExvSx_CycloBrush_",string(realV),"_T_",string(Temp),".png"];
str = join(C,"")
exportgraphics(ax ,str,'Resolution',125)

RR=(Sx*kB)';
RR(:,2)=(Ex)';
C = [pwd,"/",output,"/","ExvSx_CycloBrush_",string(realV),"_T_",string(Temp),".csv"];
str = join(C,"")
writematrix(RR,str)


plot(ExvS(:,1)*kB,ExvS(:,2),'Color',[.5 .5 .5],'LineWidth',.5)
hold on
plot(Sx*kB,Ex,'color','red','LineWidth',3)
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
