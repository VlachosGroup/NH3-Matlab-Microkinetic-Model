Surface = [];
Gas = [];
PEI = [];
TOF = [];
Temp = [];
for x=500:5:1000
Temp(end+1) = x;
[tr,sr,RR]=amm_main(x,0);
TOF_c = RR ./abyv ./(SDEN_1 + SDEN_2);
TOF(end+1,:) = reshape(TOF_c,1,3*length(TOF_c));
Surface(end+1,:)=sr{1}(end,[1:6 10]);
Gas(end+1,:)=sr{1}(end,[7:9]);
PEI(end+1,:) = RR(:,1)./(RR(:,1)+RR(:,2));
end
figure('color','w','position',[0 0 1000 1000])
subplot(2,1,1)
hold on
h=fill([0.45 0.45 0.55 0.55],[1 7 7 1],'y');
set(h,'facealpha',.1,'linestyle','none');
xlim([0,1])
set(gca,'yticklabel',{'N_2* \leftrightarrow N_2 + *';...
'2N* \leftrightarrow N_2* + *';...
'2H* \leftrightarrow H_2 + 2*';...
'NH* + * \leftrightarrow N* + H*';...
'NH_2* + * \leftrightarrow NH* + H*';...
'NH_3* + * \leftrightarrow NH_2* + H*';...
'NH_3 + * \leftrightarrow NH_3*'},'FontSize',12)
plot([0.45 0.45],[1,7],'--b')
plot([0.55 0.55],[1,7],'--b')
box on
ylabel('Reaction Step','FontSize',18)
xlabel('Partial Equilibrium Index','FontSize',18)
title({'NH_3 Decomposition on Ruthenium','Partial Equilibrium Index vs Temperature and Conversion','(\tau = 0.1 sec,  Catalyst Loading = 1500 cm^2/cm^3)'},'FontSize',18)
g=plot(PEI(1,:),[1:7],'ob', 'MarkerFacecolor','b');
g.XDataSource = 'x';
subplot(2,1,2)
plot(Temp,(1-Gas(:,3)./sum(Gas(:,1:3),2))*100)
hold on
gg=plot(Temp(1),(1-Gas(1,3)./sum(Gas(1,1:3),2))*100,'or', 'MarkerFacecolor','r');
gg.XDataSource = 'xx';
gg.YDataSource = 'yy';
xlabel('Temperature [K]','FontSize',18)
ylabel('NH_3 Conversion [%]','FontSize',18)
ylim([0 100])
set(gca, 'FontSize',12)
box on
ff = getframe(figure(1));
[im,map] = rgb2ind(ff.cdata,256,'nodither');
for z=2:size(PEI,1)
    %y=z;
    x=PEI(z,:);
    xx=Temp(z);
    yy=(1-Gas(z,3)./sum(Gas(z,1:3),2))*100;
    subplot(2,1,1)
    refreshdata
    drawnow
    ff = getframe(figure(1));
    im(:,:,1,z) = rgb2ind(ff.cdata,map,'nodither');
    %pause(0.1)
end
imwrite(im,map,'NSC_v_Temp&Conv_Ru.gif','DelayTime',0.2,'LoopCount',0);
hold off
figure(2)
BG = barh(reshape(TOF(1,:),7,3));
set(BG(3),'DisplayName','Net','FaceColor','r');
set(BG(2),'DisplayName','Reverse','FaceColor','y');
set(BG(1),'DisplayName','Forward','FaceColor','b');
set(gca,'Xscale','log', 'XMinorTick', 'off')
set(gca,'yticklabel',{'N_2* \leftrightarrow N_2 + *';...
                      '2N* \leftrightarrow N_2* + *';...
                      '2H* \leftrightarrow H_2 + 2*';...
                      'NH* + * \leftrightarrow N* + H*';...
                      'NH_2* + * \leftrightarrow NH* + H*';...
                      'NH_3* + * \leftrightarrow NH_2* + H*';...
                      'NH_3 + * \leftrightarrow NH_3*'})
xlabel('Turnover Frequency (TOF) [s^{-1}]')
ylabel('Reaction Step')
legend('Forward', 'Reverse', 'Location', 'best')
xlim([1e-30 1e10])
title(num2str(Temp(1)))
BG(1).YDataSource='y1';
BG(2).YDataSource='y2';
BG(3).YDataSource='y3';
for z=2:size(TOF,1)
    Rate=reshape(TOF(z,:),7,3);
    y1=Rate(:,1);
    y2=Rate(:,2);
    y3=Rate(:,3);
    title(num2str(Temp(z)))
    refreshdata
    drawnow
end