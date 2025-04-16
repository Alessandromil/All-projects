clc 

close all 

%STV 
[H_1, p_1] = lillietest(STV_sani.STV_30_33);
[H_2, p_2] = lillietest(STV_sani.STV_34_36); 
[H_3, p_3] = lillietest(STV_sani.STV_37_38);
H = [H_1, H_2, H_3]'; 
p = [p_1, p_2, p_3]';  
STV_sani_normalita = table(H, p);


[H_1, p_1] = lillietest(STV_malati.STV_30_33);
[H_2, p_2] = lillietest(STV_malati.STV_34_36); 
[H_3, p_3] = lillietest(STV_malati.STV_37_38);
H = [H_1, H_2, H_3]'; 
p = [p_1, p_2, p_3]'; 
STV_malati_normalita = table(H, p);
%STV ha distribuzione normale  su tutte le settimane per i sani, mentre per
%i malati la settimana 37-38 non è distribuita normalmente 

%ApEn
[H_1, p_1] = lillietest(ApEn_sani.ApEn_30_33);
[H_2, p_2] = lillietest(ApEn_sani.ApEn_34_36); 
[H_3, p_3] = lillietest(ApEn_sani.ApEn_37_38);
H = [H_1, H_2, H_3]'; 
p = [p_1, p_2, p_3]'; 
ApEn_sani_normalita = table(H, p);

[H_1, p_1] = lillietest(ApEn_malati.ApEn_30_33);
[H_2, p_2] = lillietest(ApEn_malati.ApEn_34_36); 
[H_3, p_3] = lillietest(ApEn_malati.ApEn_37_38);
H = [H_1, H_2, H_3]'; 
p = [p_1, p_2, p_3]'; 
ApEn_malati_normalita = table(H, p);
%ApEn tutti normali tranne sani (37-38)

%SampEn
[H_1, p_1] = lillietest(SampEn_sani.SampEn_30_33);
[H_2, p_2] = lillietest(SampEn_sani.SampEn_34_36); 
[H_3, p_3] = lillietest(SampEn_sani.SampEn_37_38);
H = [H_1, H_2, H_3]'; 
p = [p_1, p_2, p_3]'; 
SampEn_sani_normalita = table(H, p);

[H_1, p_1] = lillietest(SampEn_malati.SampEn_30_33);
[H_2, p_2] = lillietest(SampEn_malati.SampEn_34_36); 
[H_3, p_3] = lillietest(SampEn_malati.SampEn_37_38);
H = [H_1, H_2, H_3]'; 
p = [p_1, p_2, p_3]'; 
SampEn_malati_normalita = table(H, p);
%30-33 non normali (sani). 

%LZ2
[H_1, p_1] = lillietest(LZ2_sani.LZ2_30_33);
[H_2, p_2] = lillietest(LZ2_sani.LZ2_34_36); 
[H_3, p_3] = lillietest(LZ2_sani.LZ2_37_38);
H = [H_1, H_2, H_3]'; 
p = [p_1, p_2, p_3]'; 
LZ2_sani_normalita = table(H, p);

[H_1, p_1] = lillietest(LZ2_malati.LZ2_30_33);
[H_2, p_2] = lillietest(LZ2_malati.LZ2_34_36); 
[H_3, p_3] = lillietest(LZ2_malati.LZ2_37_38);
H = [H_1, H_2, H_3]'; 
p = [p_1, p_2, p_3]'; 
LZ2_malati_normalita = table(H, p);
%LZ2 non normali malati nella settimana 30-33. 

%LZ3
[H_1, p_1] = lillietest(LZ3_sani.LZ3_30_33);
[H_2, p_2] = lillietest(LZ3_sani.LZ3_34_36); 
[H_3, p_3] = lillietest(LZ3_sani.LZ3_37_38);
H = [H_1, H_2, H_3]'; 
p = [p_1, p_2, p_3]'; 
LZ3_sani_normalita = table(H, p);

[H_1, p_1] = lillietest(LZ3_malati.LZ3_30_33);
[H_2, p_2] = lillietest(LZ3_malati.LZ3_34_36); 
[H_3, p_3] = lillietest(LZ3_malati.LZ3_37_38);
H = [H_1, H_2, H_3]'; 
p = [p_1, p_2, p_3]'; 
LZ3_malati_normalita = table(H, p);
%LZ3 distribuito normalmente su tutte le settimane.

%HF
[H_1, p_1] = lillietest(HF_sani.HF_30_33);
[H_2, p_2] = lillietest(HF_sani.HF_34_36); 
[H_3, p_3] = lillietest(HF_sani.HF_37_38);
H = [H_1, H_2, H_3]'; 
p = [p_1, p_2, p_3]'; 
HF_sani_normalita = table(H, p);

[H_1, p_1] = lillietest(HF_malati.HF_30_33);
[H_2, p_2] = lillietest(HF_malati.HF_34_36); 
[H_3, p_3] = lillietest(HF_malati.HF_37_38);
H = [H_1, H_2, H_3]'; 
p = [p_1, p_2, p_3]'; 
HF_malati_normalita = table(H, p);
%tutti normali eccetto sani (30-33 e 37-38)

%MF
[H_1, p_1] = lillietest(MF_sani.MF_30_33);
[H_2, p_2] = lillietest(MF_sani.MF_34_36); 
[H_3, p_3] = lillietest(MF_sani.MF_37_38);
H = [H_1, H_2, H_3]'; 
p = [p_1, p_2, p_3]'; 
MF_sani_normalita = table(H, p);

[H_1, p_1] = lillietest(MF_malati.MF_30_33);
[H_2, p_2] = lillietest(MF_malati.MF_34_36); 
[H_3, p_3] = lillietest(MF_malati.MF_37_38);
H = [H_1, H_2, H_3]'; 
p = [p_1, p_2, p_3]'; 
MF_malati_normalita = table(H, p);
%MF malati e sani tutti distribuiti normalmente su tutte le settimane!

%LF
[H_1, p_1] = lillietest(LF_sani.LF_30_33);
[H_2, p_2] = lillietest(LF_sani.LF_34_36); 
[H_3, p_3] = lillietest(LF_sani.LF_37_38);
H = [H_1, H_2, H_3]'; 
p = [p_1, p_2, p_3]'; 
LF_sani_normalita = table(H, p);

[H_1, p_1] = lillietest(LF_malati.LF_30_33);
[H_2, p_2] = lillietest(LF_malati.LF_34_36); 
[H_3, p_3] = lillietest(LF_malati.LF_37_38);
H = [H_1, H_2, H_3]'; 
p = [p_1, p_2, p_3]'; 
LF_malati_normalita = table(H, p);
%LF distribuiti normalmente eccetto 30-33 sani 

%MEDIA
[H_1, p_1] = lillietest(media_sani.media_30_33);
[H_2, p_2] = lillietest(media_sani.media_34_36); 
[H_3, p_3] = lillietest(media_sani.media_37_38);
H = [H_1, H_2, H_3]'; 
p = [p_1, p_2, p_3]'; 
media_sani_normalita = table(H, p);

[H_1, p_1] = lillietest(media_malati.media_30_33);
[H_2, p_2] = lillietest(media_malati.media_34_36); 
[H_3, p_3] = lillietest(media_malati.media_37_38);
H = [H_1, H_2, H_3]'; 
p = [p_1, p_2, p_3]'; 
media_malati_normalita = table(H, p);
%media distribuita normalmente su tutti gli intervalli 

%STV ha distribuzione normale  su tutte le settimane per i sani, mentre per
%i malati la settimana 34-36 non è distribuita normalmente 
Y = [[STV_sani.STV_30_33;nan*ones(8,1)], STV_malati.STV_30_33, [STV_sani.STV_34_36;nan*ones(8,1)], STV_malati.STV_34_36, [STV_sani.STV_37_38;nan*ones(8,1)], STV_malati.STV_37_38]; 
labels = {'30-33', '30-33', '34-36', '34-36', '37-38', '37-38'};
CreaBoxplot(Y)
ylabel('STV [bpm]')

[P_30_33,H_30_33] = ranksum(STV_sani.STV_30_33, STV_malati.STV_30_33); 
[P_34_36,H_34_36] = ranksum(STV_sani.STV_34_36, STV_malati.STV_34_36); 
[P_37_38,H_37_38] = ranksum(STV_sani.STV_37_38, STV_malati.STV_37_38); 

P_STV = [P_30_33, P_34_36, P_37_38]'; 
H = [H_30_33, H_34_36, H_37_38]'; 

test_stv = table(P_STV, H);

%ApEn tutti normali tranne sani (37-38)
%Y = [[ApEn_sani.ApEn_30_33;nan*ones(8,1)], ApEn_malati.ApEn_30_33, ApEn_sani.ApEn_34_36, ApEn_malati.ApEn_34_36, ApEn_sani.ApEn_37_38, ApEn_malati.ApEn_37_38]; 
CreaBoxplot(Y)

[P_30_33,H_30_33] = ranksum(ApEn_sani.ApEn_30_33, ApEn_malati.ApEn_30_33); 
[P_34_36,H_34_36] = ranksum(ApEn_sani.ApEn_34_36, ApEn_malati.ApEn_34_36); 
[P_37_38, H_37_38] = ranksum(ApEn_sani.ApEn_37_38, ApEn_malati.ApEn_37_38); 

P_ApEn = [P_30_33, P_34_36, P_37_38]'; 
H = [H_30_33, H_34_36, H_37_38]'; 

test_ApEn = table(P_ApEn, H);

%SampEn 30-33 non normali (sani)

% Y = [SampEn_sani.SampEn_30_33, SampEn_malati.SampEn_30_33, SampEn_sani.SampEn_34_36, SampEn_malati.SampEn_34_36, SampEn_sani.SampEn_37_38, SampEn_malati.SampEn_37_38]; 
% 
% figure(3)
% h = boxplot(Y, 'labels', labels, 'Colors', 'k', 'BoxStyle','outline');
% viola = [0.4940 0.1840 0.5560]; 
% fucsia = [0.6350 0.0780 0.1840]; 
% colors = [viola; fucsia; viola; fucsia; viola; fucsia];
% m = findobj(gca,'Tag','Box');
% for j=1:length(m)
%     patch(get(m(j),'XData'),get(m(j),'YData'),colors(j, :),'FaceAlpha',.5);
% end
% set(h,{'linew'},{1})
% legend('Malati', 'Sani')
% xlabel('Settimana')
% ylabel('SampEn []')
% box on 

[P_30_33, H_30_33] = ranksum(SampEn_sani.SampEn_30_33, SampEn_malati.SampEn_30_33); 
[P_34_36,H_34_36] = ranksum(SampEn_malati.SampEn_30_33, SampEn_sani.SampEn_34_36); 
[P_37_38,H_37_38] = ranksum(SampEn_sani.SampEn_37_38, SampEn_malati.SampEn_37_38); 

P_SampEn = [P_30_33, P_34_36, P_37_38]'; 
H = [H_30_33, H_34_36, H_37_38]'; 

test_SampEn = table(P_SampEn, H);

%LZ2 tutti normali in tutte le settimane (sani e malati). 
Y = [[LZ2_sani.LZ2_30_33; nan*ones(8,1)], LZ2_malati.LZ2_30_33, [LZ2_sani.LZ2_34_36; nan*ones(8,1)], LZ2_malati.LZ2_34_36, [LZ2_sani.LZ2_37_38; nan*ones(8,1)], LZ2_malati.LZ2_37_38]; 
CreaBoxplot(Y)
ylabel('LZ2')

% figure(4)
% h = boxplot(Y, 'labels', labels, 'Colors', 'k', 'BoxStyle','outline');
% viola = [0.4940 0.1840 0.5560]; 
% fucsia = [0.6350 0.0780 0.1840]; 
% colors = [viola; fucsia; viola; fucsia; viola; fucsia];
% m = findobj(gca,'Tag','Box');
% for j=1:length(m)
%     patch(get(m(j),'XData'),get(m(j),'YData'),colors(j, :),'FaceAlpha',.5);
% end
% set(h,{'linew'},{1})
% legend('Malati', 'Sani')
% box on 
% Y = [LZ2_sani.LZ2_30_33, LZ2_malati.LZ2_30_33, LZ2_sani.LZ2_34_36, LZ2_malati.LZ2_34_36, LZ2_sani.LZ2_37_38, LZ2_malati.LZ2_37_38]; 
% labels = {'30-33', '30-33', '34-36', '34-36', '37-38', '37-38'};
% 
% figure(4)
% h = boxplot(Y, 'labels', labels, 'Colors', 'k', 'BoxStyle','outline');
% viola = [0.4940 0.1840 0.5560]; 
% fucsia = [0.6350 0.0780 0.1840]; 
% colors = [viola; fucsia; viola; fucsia; viola; fucsia];
% m = findobj(gca,'Tag','Box');
% for j=1:length(m)
%     patch(get(m(j),'XData'),get(m(j),'YData'),colors(j, :),'FaceAlpha',.5);
% end
% set(h,{'linew'},{1})
% legend('Malati', 'Sani')
% xlabel('Settimana')
% ylabel('LZ2 []')
% box on 

[P_30_33,H_30_33] = ranksum(LZ2_sani.LZ2_30_33, LZ2_malati.LZ2_30_33); 
[P_34_36,H_34_36] = ranksum(LZ2_sani.LZ2_34_36, LZ2_malati.LZ2_34_36); 
[P_37_38,H_37_38] = ranksum(LZ2_sani.LZ2_37_38, LZ2_malati.LZ2_37_38); 

P_LZ2 = [P_30_33, P_34_36, P_37_38]'; 
H = [H_30_33, H_34_36, H_37_38]'; 

test_LZ2 = table(P_LZ2, H);

%LZ3 tutti normali in tutte le settimane (sani e malati). 
Y = [[LZ3_sani.LZ3_30_33; nan*ones(8,1)], LZ3_malati.LZ3_30_33, [LZ3_sani.LZ3_34_36; nan*ones(8,1)], LZ3_malati.LZ3_34_36, [LZ3_sani.LZ3_37_38; nan*ones(8,1)], LZ3_malati.LZ3_37_38]; 
CreaBoxplot(Y)
ylabel('LZ3')
% 
% figure(5)
% h = boxplot(Y, 'labels', labels, 'Colors', 'k', 'BoxStyle','outline');
% viola = [0.4940 0.1840 0.5560]; 
% fucsia = [0.6350 0.0780 0.1840]; 
% colors = [viola; fucsia; viola; fucsia; viola; fucsia];
% m = findobj(gca,'Tag','Box');
% for j=1:length(m)
%     patch(get(m(j),'XData'),get(m(j),'YData'),colors(j, :),'FaceAlpha',.5);
% end
% set(h,{'linew'},{1})
% legend('Malati', 'Sani')
% xlabel('Settimana')
% ylabel('LZ3 []')
% box on 

[H_30_33,P_30_33] = ttest2(LZ3_sani.LZ3_30_33, LZ3_malati.LZ3_30_33); 
[H_34_36,P_34_36] = ttest2(LZ3_sani.LZ3_34_36, LZ3_malati.LZ3_34_36); 
[H_37_38,P_37_38] = ttest2(LZ3_sani.LZ3_37_38, LZ3_malati.LZ3_37_38); 

P_LZ3 = [P_30_33, P_34_36, P_37_38]'; 
H = [H_30_33, H_34_36, H_37_38]'; 

test_LZ3 = table(P_LZ3, H);

%HF tutti normali eccetto sani (30-33 e 37-38)
Y = [[HF_sani.HF_30_33;nan*ones(8,1)], HF_malati.HF_30_33, [HF_sani.HF_34_36;nan*ones(8,1)], HF_malati.HF_34_36, [HF_sani.HF_37_38;nan*ones(8,1)], HF_malati.HF_37_38]; 
CreaBoxplot(Y)
ylabel('HF [%]')
% labels = {'30-33', '30-33', '34-36', '34-36', '37-38', '37-38'};
% 
% figure(6)
% h = boxplot(Y, 'labels', labels, 'Colors', 'k', 'BoxStyle','outline');
% viola = [0.4940 0.1840 0.5560]; 
% fucsia = [0.6350 0.0780 0.1840]; 
% colors = [viola; fucsia; viola; fucsia; viola; fucsia];
% m = findobj(gca,'Tag','Box');
% for j=1:length(m)
%     patch(get(m(j),'XData'),get(m(j),'YData'),colors(j, :),'FaceAlpha',.5);
% end
% set(h,{'linew'},{1})
% xlabel('Settimane')
% ylabel('HF [%]')
% legend('Malati', 'Sani')
% box on 

[P_30_33, H_30_33] = ranksum(HF_sani.HF_30_33, HF_malati.HF_30_33); 
[P_34_36,H_34_36] = ranksum(HF_sani.HF_34_36, HF_malati.HF_34_36); 
[P_37_38, H_37_38] = ranksum(HF_sani.HF_37_38, HF_malati.HF_37_38); 

P_HF = [P_30_33, P_34_36, P_37_38]'; 
H = [H_30_33, H_34_36, H_37_38]'; 

test_HF = table(P_HF, H);

%MF malati e sani tutti distribuiti normalmente su tutte le settimane!
Y = [[MF_sani.MF_30_33; nan*ones(8,1)], MF_malati.MF_30_33, [MF_sani.MF_34_36; nan*ones(8,1)], MF_malati.MF_34_36, [MF_sani.MF_37_38; nan*ones(8,1)], MF_malati.MF_37_38]; 
CreaBoxplot(Y)
labels = {'30-33', '30-33', '34-36', '34-36', '37-38', '37-38'};
ylabel('MF')
% figure(7)
% h = boxplot(Y, 'labels', labels, 'Colors', 'k', 'BoxStyle','outline');
% viola = [0.4940 0.1840 0.5560]; 
% fucsia = [0.6350 0.0780 0.1840]; 
% colors = [viola; fucsia; viola; fucsia; viola; fucsia];
% m = findobj(gca,'Tag','Box');
% for j=1:length(m)
%     patch(get(m(j),'XData'),get(m(j),'YData'),colors(j, :),'FaceAlpha',.5);
% end
% set(h,{'linew'},{1})
% legend('Malati', 'Sani')
% xlabel('Settimana')
% ylabel('MF [%]')
% box on 

[H_30_33,P_30_33] = ttest2(MF_sani.MF_30_33, MF_malati.MF_30_33); 
[H_34_36,P_34_36] = ttest2(MF_sani.MF_34_36, MF_malati.MF_34_36); 
[H_37_38,P_37_38] = ttest2(MF_sani.MF_37_38, MF_malati.MF_37_38); 

P_MF = [P_30_33, P_34_36, P_37_38]'; 
H = [H_30_33, H_34_36, H_37_38]'; 

test_MF = table(P_MF, H);

%LF distribuiti normalmente eccetto 30-33 sani 
Y = [[LF_sani.LF_30_33; nan*ones(8,1)], LF_malati.LF_30_33, [LF_sani.LF_34_36; nan*ones(8,1)], LF_malati.LF_34_36, [LF_sani.LF_37_38; nan*ones(8,1)], LF_malati.LF_37_38]; 
CreaBoxplot(Y); 
ylabel('LF [%]')

% labels = {'30-33', '30-33', '34-36', '34-36', '37-38', '37-38'};
% 
% figure(8)
% h = boxplot(Y, 'labels', labels, 'Colors', 'k', 'BoxStyle','outline');
% viola = [0.4940 0.1840 0.5560]; 
% fucsia = [0.6350 0.0780 0.1840]; 
% colors = [viola; fucsia; viola; fucsia; viola; fucsia];
% m = findobj(gca,'Tag','Box');
% for j=1:length(m)
%     patch(get(m(j),'XData'),get(m(j),'YData'),colors(j, :),'FaceAlpha',.5);
% end
% set(h,{'linew'},{1})
% legend('Malati', 'Sani')
% xlabel('Settimane')
% ylabel('LF [%]')
% box on 

[P_30_33, H_30_33] = ranksum(LF_sani.LF_30_33, LF_malati.LF_30_33); 
[P_34_36,H_34_36] = ranksum(LF_sani.LF_34_36, LF_malati.LF_34_36); 
[P_37_38,H_37_38] = ranksum(LF_sani.LF_37_38, LF_malati.LF_37_38); 

P_LF = [P_30_33, P_34_36, P_37_38]'; 
H = [H_30_33, H_34_36, H_37_38]'; 

test_LF = table(P_LF, H);

%media distribuita normalmente su tutti gli intervalli 
Y = [[media_sani.media_30_33;nan*ones(8,1)], media_malati.media_30_33, [media_sani.media_34_36;nan*ones(8,1)], ...
    media_malati.media_34_36, [media_sani.media_37_38;nan*ones(8,1)], media_malati.media_37_38]; 
CreaBoxplot(Y)

[H_30_33,P_30_33] = ttest2(media_sani.media_30_33, media_malati.media_30_33); 
[H_34_36,P_34_36] = ttest2(media_sani.media_34_36, media_malati.media_34_36); 
[H_37_38,P_37_38] = ttest2(media_sani.media_37_38, media_malati.media_37_38); 

P_media = [P_30_33, P_34_36, P_37_38]'; 
H = [H_30_33, H_34_36, H_37_38]'; 

test_media = table(P_media, H);

close all

Y = [[ApEn_sani.ApEn_30_33; nan*ones(8,1)], ApEn_malati.ApEn_30_33, [ApEn_sani.ApEn_34_36; nan*ones(8,1)], ApEn_malati.ApEn_34_36, [ApEn_sani.ApEn_37_38; nan*ones(8,1)], ApEn_malati.ApEn_37_38]; 
CreaBoxplot(Y)
ylabel('ApEn')
