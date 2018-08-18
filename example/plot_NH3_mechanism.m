HORT = [-23.0720  -23.0503  -12.1025  -30.2070  -29.9113  -28.5829    5.3258    4.8726   -2.8737         0];
Ea = [ 0   24.5302         0   15.3072   14.1775   21.3247         0]/23.06;
e = [HORT(9) HORT(4) HORT(5)+HORT(3) HORT(6)+2*HORT(3) HORT(2)+3*HORT(3)...
    0.5*HORT(1)+3*HORT(3) 0.5*HORT(7)+3*HORT(3) 0.5*HORT(7)+1.5*HORT(8)]/23.06;
a = [Ea(7) Ea(4) Ea(5) Ea(6) Ea(2)+e(6)-e(5) Ea(1) Ea(3)];
CC = 'r';
for x=0:2:14
    plot([x x+1],[e(x/2+1)-e(1) e(x/2+1)-e(1)],CC,'LineWidth',3)
    hold on
    z1 = x+1;
    z2 = x+2;
    if x==14
        break
    end
    if a((z1+1)/2)==0
        plot([z1 z2],[e(x/2+1)-e(1) e(x/2+2)-e(1)],CC)
    else
        xx=[z1^2 z1 1;(z1+.5)^2 z1+.5 1;(z1+1)^2 z1+1 1];
        bb=[e(x/2+1) a((z1+1)/2)+e(x/2+1) e(x/2+2)]';
        coef = xx\bb;
        zz=linspace(z1,z2,100);
        plot(zz,coef(1)*zz.^2+coef(2)*zz+coef(3)-e(1),CC)
    end
end
set(gca,'xtick',0.5:1:14.5,'XTickLabelRotation',45,'xticklabel',...
    {'\bf{NH_3}','NH_3 + * \leftrightarrow NH_3*','\bf{NH_3*}',...
    'NH_3* + * \leftrightarrow NH_2* + H*','\bf{NH_2* + H*}',...
    'NH_2* + * \leftrightarrow NH* + H*','\bf{NH* + 2H*}',...
    'NH* + * \leftrightarrow N* + H*','\bf{N* + 3H*}',...
    '2N* \leftrightarrow N_2* + *','\bf{1/2N_2* + 3H*}',...
    'N_2* \leftrightarrow N_2 + *',...
    '\bf{1/2N_2 + 3H*}','2H* \leftrightarrow H_2 + 2*','\bf{1/2N_2 + 3/2H_2}'},...
    'fontsize',12)
%title({'NH_3 Decomposition Potential Energy Diagram',...
%    'Various Catalysts at Steady-State',...
%    '(T = 950[K],  P = 1[atm], \tau = 0.1[s], Loading = 1500[cm^2/cm^3], Site Density = 2.61e-9[mol/cm^2])'})
ylabel('Energy [eV]','fontsize',14)
text(3.25,a(2)+0.2+e(2),num2str(a(2),'%3.2f'),'Color',CC)
text(5.25,a(3)+0.2+e(3),num2str(a(3),'%3.2f'),'Color',CC)
text(7.25,a(4)+0.2+e(4),num2str(a(4),'%3.2f'),'Color',CC)
text(9.25,a(5)+0.2+e(5),num2str(a(5),'%3.2f'),'Color',CC)