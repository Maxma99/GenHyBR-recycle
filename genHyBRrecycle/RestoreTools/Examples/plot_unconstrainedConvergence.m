%
%  SCRIPT: plot_unconstrainedConvergence
%
%  Plot the convergence results for the test problems used in
%  IR_UnconstrainedExamples script


axes('FontSize',FS)
plot( IterInfo_ri.Enrm,'-ro','LineWidth',LW,'MarkerSize',MS), hold on
xlim([0 MaxIter+5])
plot( IterInfo_sd.Enrm,'-s','Color', myGreen,'LineWidth',LW,'MarkerSize',MS )
plot( IterInfo_lsqr.Enrm,'-bd','LineWidth',LW,'MarkerSize',MS )
plot( IterInfo_cgls.Enrm,'-m^','LineWidth',LW,'MarkerSize',MS)
plot( IterInfo_hybr.Enrm,'-kV','LineWidth',LW,'MarkerSize',MS)
legend('RI','SD', 'LSQR', 'CGLS','HyBR','Location','Best');
delete(findobj(gca,'Type','line','Marker','o'));
delete(findobj(gca,'Type','line','Marker','s'));
delete(findobj(gca,'Type','line','Marker','d'));
delete(findobj(gca,'Type','line','Marker','^'));
delete(findobj(gca,'Type','line','Marker','V'));
plot( IterInfo_ri.Enrm,'r','LineWidth',LW)
plot( IterInfo_sd.Enrm,'Color',myGreen,'LineWidth',LW)
plot( IterInfo_lsqr.Enrm,'b','LineWidth',LW )
plot( IterInfo_cgls.Enrm,'m','LineWidth',LW)
plot( IterInfo_hybr.Enrm,'k','LineWidth',LW)
plot(1:11:MaxIter,IterInfo_ri.Enrm(1:11:MaxIter),'ro','MarkerSize',MS, 'LineWidth',LW)
plot(1:11:MaxIter,IterInfo_sd.Enrm(1:11:MaxIter),'s','Color',myGreen,'MarkerSize',MS, 'LineWidth',LW)
plot(1:11:MaxIter,IterInfo_lsqr.Enrm(1:11:MaxIter),'bd','MarkerSize',MS, 'LineWidth',LW)
plot(1:11:MaxIter,IterInfo_cgls.Enrm(1:11:MaxIter),'m^','MarkerSize',MS,'LineWidth',LW)
plot(1:11:MaxIter,IterInfo_hybr.Enrm(1:11:MaxIter),'kV','MarkerSize',MS,'LineWidth',LW)
xlabel('Iterations','fontsize',FS)
ylabel('Relative error','fontsize',FS)
