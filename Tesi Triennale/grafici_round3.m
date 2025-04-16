clc 
close all
P = [P_STV'; P_media'; P_SampEn'; P_ApEn'; P_LF'; P_MF'; P_HF'; P_LZ2'; P_LZ3'];
P_prova = 1-exp(-3.*P);

P_s = zeros(size(P));

P_s(P<0.1)=1; 
P_s(P<0.05)=2;
P_s(P<0.01)=3;

figure(1)
imagesc(P_s)
colorbar
colormap("parula")
ylabels = {'STV', 'media', 'SampEn', 'ApEn', 'LF', 'MF', 'HF', 'LZ2', 'LZ3'}; 
  set(gca, 'YTick', [1; 2; 3; 4; 5; 6; 7; 8; 9], 'YTickLabel', ylabels);
  xlabels = {'', '30-33', '', '34-36', '', '37-38', ''}; 
  set(gca, 'XTick', [0.5; 1; 1.5; 2; 2.5; 3; 3.5], 'XTickLabel', xlabels);
box on 
% %stv, hf, lz3 
% 
% figure(2)
% s1 = scatter3(STV_sani.STV_30_33, MF_sani.MF_30_33, LZ3_sani.LZ3_30_33, 'filled');
% s1.MarkerFaceColor = "#0072BD"; 
% %s1.MarkerEdgeColor = "#000000";
% xlabel('STV'), ylabel('MF'), zlabel('LZ3')
% hold on 
% s2 = scatter3(STV_malati.STV_30_33, MF_malati.MF_30_33, LZ3_malati.LZ3_30_33, 'filled'); 
% s2.MarkerFaceColor = "#D95319";
% %s2.MarkerEdgeColor = "#000000";
% xlabel('STV [bpm]'), ylabel('MF [%]'), zlabel('LZ3')
% legend('Sani', 'Malati', 'Location','best')
% box on
% 
figure(2)
s1 = scatter3(STV_sani.STV_30_33, HF_sani.HF_30_33, LZ3_sani.LZ3_30_33, 'filled');
s1.MarkerFaceColor = "#0072BD"; 
%s1.MarkerEdgeColor = "#000000";
xlabel('STV'), ylabel('HF'), zlabel('LZ3')
hold on 
s2 = scatter3(STV_malati.STV_30_33, HF_malati.HF_30_33, LZ3_malati.LZ3_30_33, 'filled'); 
s2.MarkerFaceColor = "#D95319";
%s2.MarkerEdgeColor = "#000000";
xlabel('STV [bpm]'), ylabel('HF [%]'), zlabel('LZ3')
% legend('Sani', 'Malati', 'Location','best')
box on

figure(3)
s1 = scatter3(STV_sani.STV_34_36, HF_sani.HF_34_36, LZ3_sani.LZ3_34_36, 'filled');
s1.MarkerFaceColor = "#0072BD"; 
%s1.MarkerEdgeColor = "#000000";
xlabel('STV [bpm]'), ylabel('HF [%]'), zlabel('LZ3')
hold on 
s2 = scatter3(STV_malati.STV_34_36, HF_malati.HF_34_36, LZ3_malati.LZ3_34_36, 'filled'); 
s2.MarkerFaceColor = "#D95319";
%s2.MarkerEdgeColor = "#000000";
xlabel('STV [bpm]'), ylabel('HF [%]'), zlabel('LZ3')
% legend('Sani', 'Malati', 'Location','best')
box on

figure(4)
s1 = scatter3(STV_sani.STV_37_38, HF_sani.HF_37_38, LZ3_sani.LZ3_37_38, 'filled');
s1.MarkerFaceColor = "#0072BD"; 
%s1.MarkerEdgeColor = "#000000";
xlabel('STV [bpm]'), ylabel('HF [%]'), zlabel('LZ3')
hold on 
s2 = scatter3(STV_malati.STV_37_38, HF_malati.HF_37_38, LZ3_malati.LZ3_37_38, 'filled'); 
s2.MarkerFaceColor = "#D95319";
%s2.MarkerEdgeColor = "#000000";
xlabel('STV [bpm]'), ylabel('HF [%]'), zlabel('LZ3')
% legend('Sani', 'Malati', 'Location','best')
box on 

figure(5)
s1 = scatter3(STV_sani.STV_30_33, MF_sani.MF_30_33, LZ3_sani.LZ3_30_33, 'filled');
s1.MarkerFaceColor = "#0072BD"; 
%s1.MarkerEdgeColor = "#000000";
xlabel('STV'), ylabel('HF'), zlabel('LZ3')
hold on 
s2 = scatter3(STV_malati.STV_30_33, MF_malati.MF_30_33, LZ3_malati.LZ3_30_33, 'filled'); 
s2.MarkerFaceColor = "#D95319";
%s2.MarkerEdgeColor = "#000000";
xlabel('STV [bpm]'), ylabel('MF [%]'), zlabel('LZ3')
% legend('Sani', 'Malati', 'Location','best')
box on

figure(6)
s1 = scatter3(HF_sani.HF_34_36, LZ2_sani.LZ2_34_36, LZ3_sani.LZ3_34_36, 'filled');
s1.MarkerFaceColor = "#0072BD"; 
%s1.MarkerEdgeColor = "#000000";
xlabel('STV [bpm]'), ylabel('HF [%]'), zlabel('LZ3')
hold on 
s2 = scatter3(HF_malati.HF_34_36, LZ2_malati.LZ2_34_36, LZ3_malati.LZ3_34_36, 'filled'); 
s2.MarkerFaceColor = "#D95319";
%s2.MarkerEdgeColor = "#000000";
xlabel('HF [%]'), ylabel('LZ2'), zlabel('LZ3')
% legend('Sani', 'Malati', 'Location','best')
box on

figure(7)
s1 = scatter3(LF_sani.LF_37_38, HF_sani.HF_37_38, LZ3_sani.LZ3_37_38, 'filled');
s1.MarkerFaceColor = "#0072BD"; 
%s1.MarkerEdgeColor = "#000000";
xlabel('STV [bpm]'), ylabel('HF [%]'), zlabel('LZ3')
hold on 
s2 = scatter3(LF_malati.LF_37_38, HF_malati.HF_37_38, LZ3_malati.LZ3_37_38, 'filled'); 
s2.MarkerFaceColor = "#D95319";
%s2.MarkerEdgeColor = "#000000";
xlabel('LF [%]'), ylabel('HF [%]'), zlabel('LZ3')
% legend('Sani', 'Malati', 'Location','best')
box on 

% 
% 
% %STV ha distribuzione normale  su tutte le settimane per i sani, mentre per
% %i malati la settimana 34-36 non è distribuita normalmente 
% mean_1 = mean(STV_sani.STV_30_33); 
% dev = std(STV_sani.STV_30_33); 
% stv_par_s1 = [mean_1; dev]; 
% 
% mean_1 = mean(STV_malati.STV_30_33); 
% dev = std(STV_malati.STV_30_33); 
% stv_par_m1 = [mean_1; dev]; 
% 
% mean_1 = mean(STV_sani.STV_34_36); 
% dev = std(STV_sani.STV_34_36); 
% stv_par_s2 = [mean_1; dev]; 
% 
% dev = quantile(STV_malati.STV_34_36, [0.25, 0.5, 0.75]); 
% stv_par_m2 = dev; 
% 
% mean_1 = mean(STV_sani.STV_37_38); 
% dev = std(STV_sani.STV_37_38); 
% stv_par_s3 = [mean_1; dev]; 
% 
% mean_1 = mean(STV_malati.STV_37_38); 
% dev = std(STV_malati.STV_37_38); 
% stv_par_m3 = [mean_1; dev]; 
% 
% %ApEn tutti normali tranne sani (37-38)
% mean_1 = mean(ApEn_sani.ApEn_30_33); 
% dev = std(ApEn_sani.ApEn_30_33); 
% ApEn_par_s1 = [mean_1; dev];
% 
% mean_1 = mean(ApEn_sani.ApEn_34_36); 
% dev = std(ApEn_sani.ApEn_34_36); 
% ApEn_par_s2 = [mean_1; dev];
% 
% dev = quantile(ApEn_sani.ApEn_37_38, [0.25; 0.5; 0.75]); 
% ApEn_par_s3 = dev;
% 
% mean_1 = mean(ApEn_malati.ApEn_30_33); 
% dev = std(ApEn_malati.ApEn_30_33); 
% ApEn_par_m1 = [mean_1; dev];
% 
% mean_1 = mean(ApEn_malati.ApEn_34_36); 
% dev = std(ApEn_malati.ApEn_34_36); 
% ApEn_par_m2 = [mean_1; dev];
% 
% mean_1 = mean(ApEn_malati.ApEn_37_38); 
% dev = std(ApEn_malati.ApEn_37_38); 
% ApEn_par_m3 = [mean_1; dev];
% 
% %SampEn 30-33 non normali (sani). 
% dev = quantile(SampEn_sani.SampEn_30_33, [0.25; 0.5; 0.75]); 
% SampEn_par_s1 = dev;
% 
% mean_1 = mean(SampEn_sani.SampEn_34_36); 
% dev = std(SampEn_sani.SampEn_34_36); 
% SampEn_par_s2 = [mean_1; dev];
% 
% mean_1 = mean(SampEn_sani.SampEn_37_38); 
% dev = std(SampEn_sani.SampEn_37_38); 
% SampEn_par_s3 = [mean_1; dev];
% 
% mean_1 = mean(SampEn_malati.SampEn_30_33); 
% dev = std(SampEn_malati.SampEn_30_33); 
% SampEn_par_m1 = [mean_1; dev];
% 
% mean_1 = mean(SampEn_malati.SampEn_34_36); 
% dev = std(SampEn_malati.SampEn_34_36); 
% SampEn_par_m2 = [mean_1; dev]; 
% 
% mean_1 = mean(SampEn_malati.SampEn_37_38); 
% dev = std(SampEn_malati.SampEn_37_38); 
% SampEn_par_m3 = [mean_1; dev];
% 
% %LZ2 tutti normali in tutte le settimane (sani e malati). 
% mean_1 = mean(LZ2_sani.LZ2_30_33);
% dev = std(LZ2_sani.LZ2_30_33); 
% LZ2_par_s1 = [mean_1; dev];
% 
% mean_1 = mean(LZ2_sani.LZ2_34_36); 
% dev = std(LZ2_sani.LZ2_34_36); 
% LZ2_par_s2 = [mean_1; dev];
% 
% mean_1 = mean(LZ2_sani.LZ2_37_38); 
% dev = std(LZ2_sani.LZ2_37_38); 
% LZ2_par_s3 = [mean_1; dev];
% 
% mean_1 = mean(LZ2_malati.LZ2_30_33); 
% dev = std(LZ2_malati.LZ2_30_33); 
% LZ2_par_m1 = [mean_1; dev];
% 
% mean_1 = mean(LZ2_malati.LZ2_34_36); 
% dev = std(LZ2_malati.LZ2_34_36); 
% LZ2_par_m2 = [mean_1; dev]; 
% 
% mean_1 = mean(LZ2_malati.LZ2_37_38); 
% dev = std(LZ2_malati.LZ2_37_38); 
% LZ2_par_m3 = [mean_1; dev];
% 
% %LZ3 distribuito normalmente su tutte le settimane a parte 34-36 malati 
% mean_1 = mean(LZ3_sani.LZ3_30_33);
% dev = std(LZ3_sani.LZ3_30_33); 
% LZ3_par_s1 = [mean_1; dev];
% 
% mean_1 = mean(LZ3_sani.LZ3_34_36); 
% dev = std(LZ3_sani.LZ3_34_36); 
% LZ3_par_s2 = [mean_1; dev];
% 
% mean_1 = mean(LZ3_sani.LZ3_37_38); 
% dev = std(LZ3_sani.LZ3_37_38); 
% LZ3_par_s3 = [mean_1; dev];
%  
% dev = quantile(LZ3_malati.LZ3_30_33, [0.25; 0.5; 0.75]); 
% LZ3_par_m1 = dev;
% 
% mean_1 = mean(LZ3_malati.LZ3_34_36); 
% dev = std(LZ3_malati.LZ3_34_36); 
% LZ3_par_m2 = [mean_1; dev]; 
% 
% mean_1 = mean(LZ3_malati.LZ3_37_38); 
% dev = std(LZ3_malati.LZ3_37_38); 
% LZ3_par_m3 = [mean_1; dev];
% 
% %HF tutti normali eccetto sani (30-33 e 37-38)
% 
% dev = quantile(HF_sani.HF_30_33, [0.25; 0.5; 0.75]); 
% HF_par_s1 = dev;
% 
% mean_1 = mean(HF_sani.HF_34_36); 
% dev = std(HF_sani.HF_34_36); 
% HF_par_s2 = [mean_1; dev];
%  
% dev = quantile(HF_sani.HF_37_38, [0.25; 0.5; 0.75]); 
% HF_par_s3 = dev;
% 
% mean_1 = mean(HF_malati.HF_30_33);
% dev = std(HF_malati.HF_30_33); 
% HF_par_m1 = [mean_1; dev];
% 
% mean_1 = mean(HF_malati.HF_34_36); 
% dev = std(HF_malati.HF_34_36); 
% HF_par_m2 = [mean_1; dev]; 
% 
% mean_1 = mean(HF_malati.HF_37_38); 
% dev = std(HF_malati.HF_37_38); 
% HF_par_m3 = [mean_1; dev];
% 
% %MF malati e sani tutti distribuiti normalmente su tutte le settimane!
% mean_1 = mean(MF_sani.MF_30_33);
% dev = std(MF_sani.MF_30_33); 
% MF_par_s1 = [mean_1; dev];
% 
% mean_1 = mean(MF_sani.MF_34_36); 
% dev = std(MF_sani.MF_34_36); 
% MF_par_s2 = [mean_1; dev];
% 
% mean_1 = mean(MF_sani.MF_37_38); 
% dev = std(MF_sani.MF_37_38); 
% MF_par_s3 = [mean_1; dev];
% 
% mean_1 = mean(MF_malati.MF_30_33);
% dev = std(MF_malati.MF_30_33); 
% MF_par_m1 = [mean_1; dev];
% 
% mean_1 = mean(MF_malati.MF_34_36); 
% dev = std(MF_malati.MF_34_36); 
% MF_par_m2 = [mean_1; dev]; 
% 
% mean_1 = mean(MF_malati.MF_37_38); 
% dev = std(MF_malati.MF_37_38); 
% MF_par_m3 = [mean_1; dev];
% 
% %LF distribuiti normalmente eccetto 30-33 sani 
% dev = quantile(LF_sani.LF_30_33, [0.25; 0.5; 0.75]); 
% LF_par_s1 = dev;
% 
% mean_1 = mean(LF_sani.LF_34_36); 
% dev = std(LF_sani.LF_34_36); 
% LF_par_s2 = [mean_1; dev];
% 
% mean_1 = mean(LF_sani.LF_37_38); 
% dev = std(LF_sani.LF_37_38); 
% LF_par_s3 = [mean_1; dev];
% 
% mean_1 = mean(LF_malati.LF_30_33);
% dev = std(LF_malati.LF_30_33); 
% LF_par_m1 = [mean_1; dev];
% 
% mean_1 = mean(LF_malati.LF_34_36); 
% dev = std(LF_malati.LF_34_36); 
% LF_par_m2 = [mean_1; dev]; 
% 
% mean_1 = mean(LF_malati.LF_37_38); 
% dev = std(LF_malati.LF_37_38); 
% LF_par_m3 = [mean_1; dev];
% 
% %media distribuita normalmente su tutti gli intervalli 
% mean_1 = mean(media_sani.media_30_33);
% dev = std(media_sani.media_30_33); 
% media_par_s1 = [mean_1; dev];
% 
% mean_1 = mean(media_sani.media_34_36); 
% dev = std(media_sani.media_34_36); 
% media_par_s2 = [mean_1; dev];
% 
% mean_1 = mean(media_sani.media_37_38); 
% dev = std(media_sani.media_37_38); 
% media_par_s3 = [mean_1; dev];
% 
% mean_1 = mean(media_malati.media_30_33);
% dev = std(media_malati.media_30_33); 
% media_par_m1 = [mean_1; dev];
% 
% mean_1 = mean(media_malati.media_34_36); 
% dev = std(media_malati.media_34_36); 
% media_par_m2 = [mean_1; dev]; 
% 
% mean_1 = mean(media_malati.media_37_38); 
% dev = std(media_malati.media_37_38); 
% media_par_m3 = [mean_1; dev];

ApEn_FrS = [ApEn_sani.ApEn_30_33, ApEn_sani.ApEn_34_36, ApEn_sani.ApEn_37_38];
ApEn_FrM = [ApEn_malati.ApEn_30_33, ApEn_malati.ApEn_34_36, ApEn_malati.ApEn_37_38]; 

ApEn_Fr = [ApEn_FrS;ApEn_FrM];

SampEn_Frs = [SampEn_sani.SampEn_30_33, SampEn_sani.SampEn_34_36, SampEn_sani.SampEn_37_38]; 
SampEn_Frm = [SampEn_malati.SampEn_30_33, SampEn_malati.SampEn_34_36, SampEn_malati.SampEn_37_38]; 
HF_Frs = [HF_sani.HF_30_33, HF_sani.HF_34_36, HF_sani.HF_37_38]; 
HF_Frm = [HF_malati.HF_30_33, HF_malati.HF_34_36, HF_malati.HF_37_38];
MF_Frs = [MF_sani.MF_30_33, MF_sani.MF_34_36, MF_sani.MF_37_38]; 
MF_Frm = [MF_malati.MF_30_33, MF_malati.MF_34_36, MF_malati.MF_37_38]; 
LF_Frs = [LF_sani.LF_30_33, LF_sani.LF_34_36, LF_sani.LF_37_38]; 
LF_Frm = [LF_malati.LF_30_33, LF_malati.LF_34_36, LF_malati.LF_37_38]; 
STV_Frs = [STV_sani.STV_30_33, STV_sani.STV_34_36, STV_sani.STV_37_38]; 
STV_Frm = [STV_malati.STV_30_33, STV_malati.STV_34_36, STV_malati.STV_37_38]; 
LZ2_Frs = [LZ2_sani.LZ2_30_33, LZ2_sani.LZ2_34_36, LZ2_sani.LZ2_37_38];
LZ2_Frm = [LZ2_malati.LZ2_30_33, LZ2_malati.LZ2_34_36, LZ2_malati.LZ2_37_38]; 
LZ3_Frs = [LZ3_sani.LZ3_30_33, LZ3_sani.LZ3_34_36, LZ3_sani.LZ3_37_38]; 
LZ3_Frm = [LZ3_malati.LZ3_30_33, LZ3_malati.LZ3_34_36, LZ3_malati.LZ3_37_38];
media_FrS = [media_sani.media_30_33, media_sani.media_34_36, media_sani.media_37_38]; 
media_FrM = [media_malati.media_30_33, media_malati.media_34_36, media_malati.media_37_38];

[P, table, statsS] = friedman(ApEn_FrS); 
[P, table, statsM] = friedman(ApEn_FrM);
% figure(1)
% multcompare(stats)
% title('ApEn')

[P, table, statsS] = friedman(STV_Frs);
[P, table, statsM] = friedman(STV_Frm);
% figure(2)
% multcompare(stats)
% title('STV')

[P, table, statsS] = friedman(SampEn_Frs);
[P, table, statsM] = friedman(SampEn_Frm);

[P, table, statsS] = friedman(HF_Frs);
[P, table, statsM] = friedman(HF_Frm);

[P, table, statsS] = friedman(MF_Frs); 
[P, table, statsM] = friedman(MF_Frm); 

[P, table, statsS] = friedman(LF_Frs); 
[P, table, statsM] = friedman(LF_Frm); 

[P, table, statsS] = friedman(LZ2_Frs); 
[P, table, statsM] = friedman(LZ2_Frm); 

[P, table, statsS] = friedman(LZ3_Frs); 
[P, table, statsM] = friedman(LZ3_Frm); 

[P, table, statsS] = friedman(media_FrS);
[P, table, statsM] = friedman(media_FrM);

figure
boxplot(media_FrS)
hold on
plot(media_FrS')

figure(1)
plot(ApEn_FrM')


