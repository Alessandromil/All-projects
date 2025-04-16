clc
clear all
close all
load('data_Team_6.mat')


%% Graficos SV

figure(1) 
%SV_16
subplot(3,1,1) 
plot(t_16,Q_16) 
xlabel('t') 
ylabel('Q') 
axis([0 1 -50 200]) 
title('SV_16') 
%SV_17
subplot(3,1,2) 
plot(t_17,Q_17) 
xlabel('t') 
ylabel('Q') 
axis([0 1 -50 500]) 
title('SV_17') 
%SV_18
subplot(3,1,3) 
plot(t_18,Q_18) 
xlabel('t') 
ylabel('Q') 
axis([0 1 -50 500]) 
title('SV_18') 

%% Calcula SV

integralSV16 = trapz(t_16,Q_16); 
SV_16 = (integralSV16)/1000
HR_16 = 60/1.024
CO_16 = HR_16 * SV_16
Pa_16 = mean(P_16);
R_16 = Pa_16 / CO_16


integralSV17 = trapz(t_17,Q_17); 
SV_17 = (integralSV17)/1000
HR_17 = 60/0.822
CO_17 = HR_17 * SV_17
Pa_17 = mean(P_17);
R_17 = Pa_17 / CO_17


integralSV18 = trapz(t_18,Q_18); 
SV_18 = (integralSV18)/1000
HR_18 = 60/0.97
CO_18 = HR_18 * SV_18
Pa_18 = mean(P_18);
R_18 = Pa_18 / CO_18
%16 is a child, 17 and 18 are adults

%% Pin_16

cte_16=40-R1_16*Q_16(1)-Pout_16;
QQ_16 =[Q_16;Q_16;Q_16;Q_16;Q_16;Q_16;Q_16;Q_16];
PP_16 =[P_16;P_16;P_16;P_16;P_16;P_16;P_16;P_16];

tt_16 =zeros(length(QQ_16));
for k = 2:length(tt_16)
    tt_16(k)=tt_16(k-1)+0.001;
end
Pin_16=zeros(length(tt_16),1);
for i=1:length(tt_16)

    Pin_16(i)=cte_16*(exp(-(tt_16(i)-tt_16(1))/(R2_16*C_16)))+R1_16*QQ_16(i)+Pout_16;
    tiempo=zeros(i,1);

    caudal=tiempo;

    for j=1:i

        tiempo(j,1)=tt_16(j,1);
        caudal(j,1)=exp(tiempo(j,1)/(R2_16*C_16))*QQ_16(j,1);

    end

    if i==1

        Pin_16(i,1)=Pin_16(i,1);

    else

        Pin_16(i,1)=Pin_16(i,1)+(exp(-tt_16(i,1)/(R2_16*C_16))/C_16)*trapz(tiempo,caudal);

    end

end

Pin_16ult=zeros(length(t_16),1);
l=1;
for i=length(tt_16)-length(t_16)+1:1:length(tt_16)

    Pin_16ult(l)=Pin_16(i);
    l=l+1;
end


% figure(2)
% subplot(3,1,1)
% plot(tt_16,Pin_16,'b')
% hold on;
% plot(tt_16,PP_16,'r')
% hold off;
% title('Pin 16')

%Pin_17

cte_17=80-R1_17*Q_17(1)-Pout_17;
QQ_17 = [Q_17; Q_17; Q_17; Q_17; Q_17; Q_17; Q_17; Q_17];
PP_17 = [P_17; P_17; P_17; P_17; P_17; P_17; P_17; P_17];

tt_17 = zeros(length(QQ_17));
for k = 2:length(tt_17)
    tt_17(k) = tt_17(k-1) + 0.001;
end

Pin_17 = zeros(length(tt_17), 1);
for i = 1:length(tt_17)

    Pin_17(i) = cte_17 * (exp(-(tt_17(i) - tt_17(1)) / (R2_17 * C_17))) + R1_17 * QQ_17(i) + Pout_17;
    tiempo2 = zeros(i, 1);

    caudal2 = tiempo2;

    for j = 1:i

        tiempo2(j, 1) = tt_17(j, 1);
        caudal2(j, 1) = exp(tiempo2(j, 1) / (R2_17 * C_17)) * QQ_17(j, 1);

    end

    if i == 1
        Pin_17(i, 1) = Pin_17(i, 1);
    else
        Pin_17(i, 1) = Pin_17(i, 1) + (exp(-tt_17(i, 1) / (R2_17 * C_17)) / C_17) * trapz(tiempo2, caudal2);
    end

end


Pin_17ult=zeros(length(t_17),1);
l=1;
for i=length(tt_17)-length(t_17)+1:1:length(tt_17)

    Pin_17ult(l)=Pin_17(i);
    l=l+1;
end


% subplot(3, 1, 2)
% plot(tt_17,Pin_17,'b')
% hold on;
% plot(tt_17,PP_17,'r')
% hold off;
% title('Pin 17')

%Pin_18
cte_18=80-R1_18*Q_18(1)-Pout_18;
QQ_18 = [Q_18; Q_18; Q_18; Q_18; Q_18; Q_18; Q_18; Q_18];
PP_18 = [P_18; P_18; P_18; P_18; P_18; P_18; P_18; P_18];

tt_18 = zeros(length(QQ_18));

for k = 2:length(tt_18)
    tt_18(k) = tt_18(k-1) + 0.001;
end

Pin_18 = zeros(length(tt_18), 1);
for i = 1:length(tt_18)

    Pin_18(i) = cte_18 * (exp(-(tt_18(i) - tt_18(1)) / (R2_18 * C_18))) + R1_18 * QQ_18(i) + Pout_18;
    tiempo3 = zeros(i, 1);

    caudal3 = tiempo3;

    for j = 1:i

        tiempo3(j, 1) = tt_18(j, 1);
        caudal3(j, 1) = exp(tiempo3(j, 1) / (R2_18 * C_18)) * QQ_18(j, 1);

    end

    if i == 1
        Pin_18(i, 1) = Pin_18(i, 1);
    else
        Pin_18(i, 1) = Pin_18(i, 1) + (exp(-tt_18(i, 1) / (R2_18 * C_18)) / C_18) * trapz(tiempo3, caudal3);
    end

end

Pin_18ult=zeros(length(t_18),1);
l=1;
for i=length(tt_18)-length(t_18)+1:1:length(tt_18)

    Pin_18ult(l)=Pin_18(i);
    l=l+1;
end

% subplot(3,1,3)
% plot(tt_18,Pin_18,'b')
% hold on;
% plot(tt_18,PP_18,'r')
% hold off;
% title('Pin 18')


figure(3)
subplot(3,1,1)
plot(t_16,Pin_16ult,'b')
hold on;
plot(t_16,P_16,'r')
hold off;
title('Pin 16')

subplot(3, 1, 2)
plot(t_17, Pin_17ult, 'b')
hold on;
plot(t_17, P_17, 'r')
hold off;
title('Pin 17')

subplot(3, 1, 3)
plot(t_18, Pin_18ult, 'b')
hold on;
plot(t_18, P_18, 'r')
hold off;
title('Pin 18')

%% Errors

for i=1:length(t_16)

Errorepp_16=(1/length(t_16))*sqrt(sum(((Pin_16ult(i)-P_16(i))/P_16(i)).^2))*100;

end

Errorepp_16 = Errorepp_16
Erroreavg_16=abs(((mean(Pin_16ult)-mean(P_16))/mean(P_16))*100)
Erooresys_16=abs(((max(Pin_16ult)-max(P_16))/max(P_16))*100)
Erroredia_16=abs(((min(Pin_16ult)-min(P_16))/min(P_16))*100)




for i=1:length(t_17)

Errorepp_17=(1/length(t_17))*sqrt(sum(((Pin_17ult(i)-P_17(i))/P_17(i)).^2))*100;

end

Errorepp_17 = Errorepp_17
Erroreavg_17=abs(((mean(Pin_17)-mean(P_17))/mean(P_17))*100)
Erooresys_17=abs(((max(Pin_17ult)-max(P_17))/max(P_17))*100)
Erroredia_17=abs(((min(Pin_17ult)-min(P_17))/min(P_17))*100)



for i=1:length(t_18)

Errorepp_18=(1/length(t_18))*sqrt(sum(((Pin_18ult(i)-P_18(i))/P_18(i)).^2))*100;

end

Errorepp_18 = Errorepp_18
Erroreavg_18=abs(((mean(Pin_18ult)-mean(P_18))/mean(P_18))*100)
Erooresys_18=abs(((max(Pin_18ult)-max(P_18))/max(P_18))*100)
Erroredia_18=abs(((min(Pin_18ult)-min(P_18))/min(P_18))*100)

%% Point 3

R2_16inc=1.25*R2_16;
R2_16dec=0.75*R2_16;

cte_16dec=40-R1_16*Q_16(1)-Pout_16;
cte_16inc=40-R1_16*Q_16(1)-Pout_16;
QQ_16 =[Q_16;Q_16;Q_16;Q_16;Q_16;Q_16;Q_16;Q_16];
PP_16 =[P_16;P_16;P_16;P_16;P_16;P_16;P_16;P_16];

tt_16 =zeros(length(QQ_16));
for k = 2:length(tt_16)
    tt_16(k)=tt_16(k-1)+0.001;
end
Pin_16inc=zeros(length(tt_16),1);
for i=1:length(tt_16)

    Pin_16inc(i)=cte_16inc*(exp(-(tt_16(i)-tt_16(1))/(R2_16inc*C_16)))+R1_16*QQ_16(i)+Pout_16;
    tiempo=zeros(i,1);

    caudal=tiempo;

    for j=1:i

        tiempo(j,1)=tt_16(j,1);
        caudal(j,1)=exp(tiempo(j,1)/(R2_16inc*C_16))*QQ_16(j,1);

    end

    if i==1

        Pin_16inc(i,1)=Pin_16inc(i,1);

    else

        Pin_16inc(i,1)=Pin_16inc(i,1)+(exp(-tt_16(i,1)/(R2_16inc*C_16))/C_16)*trapz(tiempo,caudal);

    end

end

figure(4)
subplot(2,1,1)
plot(tt_16,Pin_16inc,'b')
hold on;
plot(tt_16,PP_16,'r')
hold off;
title('Pin 16 inc')

tt_16 =zeros(length(QQ_16));
for k = 2:length(tt_16)
    tt_16(k)=tt_16(k-1)+0.001;
end
Pin_16dec=zeros(length(tt_16),1);
for i=1:length(tt_16)

    Pin_16dec(i)=cte_16dec*(exp(-(tt_16(i)-tt_16(1))/(R2_16dec*C_16)))+R1_16*QQ_16(i)+Pout_16;
    tiempo=zeros(i,1);

    caudal=tiempo;

    for j=1:i

        tiempo(j,1)=tt_16(j,1);
        caudal(j,1)=exp(tiempo(j,1)/(R2_16dec*C_16))*QQ_16(j,1);

    end

    if i==1

        Pin_16dec(i,1)=Pin_16dec(i,1);

    else

        Pin_16dec(i,1)=Pin_16dec(i,1)+(exp(-tt_16(i,1)/(R2_16dec*C_16))/C_16)*trapz(tiempo,caudal);

    end

end

subplot(2,1,2)
plot(tt_16,Pin_16dec,'b')
hold on;
plot(tt_16,PP_16,'r')
hold off;
title('Pin 16 dec')
%% 17

R2_17inc=1.25*R2_17;
R2_17dec=0.75*R2_17;

cte_17inc=100-R1_17*Q_17(1)-Pout_17;
cte_17dec=70-R1_17*Q_17(1)-Pout_17;
QQ_17 =[Q_17;Q_17;Q_17;Q_17;Q_17;Q_17;Q_17;Q_17];
PP_17 =[P_17;P_17;P_17;P_17;P_17;P_17;P_17;P_17];

tt_17 =zeros(length(QQ_17));
for k = 2:length(tt_17)
    tt_17(k)=tt_17(k-1)+0.001;
end
Pin_17inc=zeros(length(tt_17),1);
for i=1:length(tt_17)

    Pin_17inc(i)=cte_17inc*(exp(-(tt_17(i)-tt_17(1))/(R2_17inc*C_17)))+R1_17*QQ_17(i)+Pout_17;
    tiempo=zeros(i,1);

    caudal=tiempo;

    for j=1:i

        tiempo(j,1)=tt_17(j,1);
        caudal(j,1)=exp(tiempo(j,1)/(R2_17inc*C_17))*QQ_17(j,1);

    end

    if i==1

        Pin_17inc(i,1)=Pin_17inc(i,1);

    else

        Pin_17inc(i,1)=Pin_17inc(i,1)+(exp(-tt_17(i,1)/(R2_17inc*C_17))/C_17)*trapz(tiempo,caudal);

    end

end

figure(5)
subplot(2,1,1)
plot(tt_17,Pin_17inc,'b')
hold on;
plot(tt_17,PP_17,'r')
hold off;
title('Pin 17 inc')

tt_17 =zeros(length(QQ_17));
for k = 2:length(tt_17)
    tt_17(k)=tt_17(k-1)+0.001;
end
Pin_17dec=zeros(length(tt_17),1);
for i=1:length(tt_17)

    Pin_17dec(i)=cte_17dec*(exp(-(tt_17(i)-tt_17(1))/(R2_17dec*C_17)))+R1_17*QQ_17(i)+Pout_17;
    tiempo=zeros(i,1);

    caudal=tiempo;

    for j=1:i

        tiempo(j,1)=tt_17(j,1);
        caudal(j,1)=exp(tiempo(j,1)/(R2_17dec*C_17))*QQ_17(j,1);

    end

    if i==1

        Pin_17dec(i,1)=Pin_17dec(i,1);

    else

        Pin_17dec(i,1)=Pin_17dec(i,1)+(exp(-tt_17(i,1)/(R2_17dec*C_17))/C_17)*trapz(tiempo,caudal);

    end

end

subplot(2,1,2)
plot(tt_17,Pin_17dec,'b')
hold on;
plot(tt_17,PP_17,'r')
hold off;
title('Pin 17 dec')

%% 18
R2_18inc=1.25*R2_18;
R2_18dec=0.75*R2_18;

cte_18inc=105-R1_18*Q_18(1)-Pout_18;
cte_18dec=73-R1_18*Q_18(1)-Pout_18;
QQ_18 =[Q_18;Q_18;Q_18;Q_18;Q_18;Q_18;Q_18;Q_18];
PP_18 =[P_18;P_18;P_18;P_18;P_18;P_18;P_18;P_18];

tt_18 =zeros(length(QQ_18));
for k = 2:length(tt_18)
    tt_18(k)=tt_18(k-1)+0.001;
end
Pin_18inc=zeros(length(tt_18),1);
for i=1:length(tt_18)

    Pin_18inc(i)=cte_18inc*(exp(-(tt_18(i)-tt_18(1))/(R2_18inc*C_18)))+R1_18*QQ_18(i)+Pout_18;
    tiempo=zeros(i,1);

    caudal=tiempo;

    for j=1:i

        tiempo(j,1)=tt_18(j,1);
        caudal(j,1)=exp(tiempo(j,1)/(R2_18inc*C_18))*QQ_18(j,1);

    end

    if i==1

        Pin_18inc(i,1)=Pin_18inc(i,1);

    else

        Pin_18inc(i,1)=Pin_18inc(i,1)+(exp(-tt_18(i,1)/(R2_18inc*C_18))/C_18)*trapz(tiempo,caudal);

    end

end

figure(6)
subplot(2,1,1)
plot(tt_18,Pin_18inc,'b')
hold on;
plot(tt_18,PP_18,'r')
hold off;
title('Pin 18 inc')

tt_18 =zeros(length(QQ_18));
for k = 2:length(tt_18)
    tt_18(k)=tt_18(k-1)+0.001;
end
Pin_18dec=zeros(length(tt_18),1);
for i=1:length(tt_18)

    Pin_18dec(i)=cte_18dec*(exp(-(tt_18(i)-tt_18(1))/(R2_18dec*C_18)))+R1_18*QQ_18(i)+Pout_18;
    tiempo=zeros(i,1);

    caudal=tiempo;

    for j=1:i

        tiempo(j,1)=tt_18(j,1);
        caudal(j,1)=exp(tiempo(j,1)/(R2_18dec*C_18))*QQ_18(j,1);

    end

    if i==1

        Pin_18dec(i,1)=Pin_18dec(i,1);

    else

        Pin_18dec(i,1)=Pin_18dec(i,1)+(exp(-tt_18(i,1)/(R2_18dec*C_18))/C_18)*trapz(tiempo,caudal);

    end

end

subplot(2,1,2)
plot(tt_18,Pin_18dec,'b')
hold on;
plot(tt_18,PP_18,'r')
hold off;
title('Pin 18 dec')

%% Ult

Pin_16dec_ult=zeros(length(t_16),1);
l=1;
for i=length(tt_16)-length(t_16)+1:1:length(tt_16)

    Pin_16dec_ult(l)=Pin_16dec(i);
    l=l+1;
end
Pin_16inc_ult=zeros(length(t_16),1);
l=1;
for i=length(tt_16)-length(t_16)+1:1:length(tt_16)

    Pin_16inc_ult(l)=Pin_16inc(i);
    l=l+1;
end




Pin_17dec_ult=zeros(length(t_17),1);
l=1;
for i=length(tt_17)-length(t_17)+1:1:length(tt_17)

    Pin_17dec_ult(l)=Pin_17dec(i);
    l=l+1;
end
Pin_17inc_ult=zeros(length(t_17),1);
l=1;
for i=length(tt_17)-length(t_17)+1:1:length(tt_17)

    Pin_17inc_ult(l)=Pin_17inc(i);
    l=l+1;
end






Pin_18dec_ult=zeros(length(t_18),1);
l=1;
for i=length(tt_18)-length(t_18)+1:1:length(tt_18)

    Pin_18dec_ult(l)=Pin_18dec(i);
    l=l+1;
end
Pin_18inc_ult=zeros(length(t_18),1);
l=1;
for i=length(tt_18)-length(t_18)+1:1:length(tt_18)

    Pin_18inc_ult(l)=Pin_18inc(i);
    l=l+1;
end






figure(7)
subplot(3,1,1)
plot(t_16,Pin_16dec_ult,'b')
hold on;
plot(t_16,P_16,'r')
hold off;
title('Pin 16 decreased')

subplot(3, 1, 2)
plot(t_17, Pin_17dec_ult, 'b')
hold on;
plot(t_17, P_17, 'r')
hold off;
title('Pin 17 decreased')

subplot(3, 1, 3)
plot(t_18, Pin_18dec_ult, 'b')
hold on;
plot(t_18, P_18, 'r')
hold off;
title('Pin 18 decreased')



figure(8)
subplot(3,1,1)
plot(t_16,Pin_16inc_ult,'b')
hold on;
plot(t_16,P_16,'r')
hold off;
title('Pin 16 increased')

subplot(3, 1, 2)
plot(t_17, Pin_17inc_ult, 'b')
hold on;
plot(t_17, P_17, 'r')
hold off;
title('Pin 17 increased')

subplot(3, 1, 3)
plot(t_18, Pin_18inc_ult, 'b')
hold on;
plot(t_18, P_18, 'r')
hold off;
title('Pin 18 increased')
