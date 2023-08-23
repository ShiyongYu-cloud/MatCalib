function [] = plot_difference(diff_pdfs)
%% function for plotting the pdf and hpd regions of the difference between two modeled ages
%INPUT
%diff_pdfs: matrix containing the estimated pdf and its 68.2% and 95.4%
%higher probabilities
%%
X = diff_pdfs(:,1);
Xprob = diff_pdfs(:,2);
X68p2 = diff_pdfs(:,3);
X95p4 = diff_pdfs(:,4);
X0 = X(Xprob>0.000001);
P0 = Xprob(Xprob>0.000001);
X1 = X(X95p4>0.000001);
P1 = X95p4(X95p4>0.000001);
X2 = X(X68p2>0.000001);
P2 = X68p2(X68p2>0.000001);
%% plot pdf
fill([X0(1);X0;X0(end);X0(1)],[0;P0;0;0],[0.301 0.745 0.933]);
hold on
%% plot the 95.4% hpd region
fill([X1(1);X1;X1(end);X1(1)],[0;P1;0;0],[0.929 0.694 0.125]);
hold on
%% plot the 68.2% hpd region
fill([X2(1);X2;X2(end);X2(1)],[0;P2;0;0],[0.635 0.078 0.184]);
grid on;
xlabel('Span (years)');
ylabel('Probability');
legend('Posterior probability','95.4 % HPD regions','68.2 % HPD regions'); 
set(gca,'XMinorTick','on');
set(gca,'TickDir','out');
end