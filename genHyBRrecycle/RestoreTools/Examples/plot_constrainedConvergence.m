%
%  SCRIPT: plot_constrainedConvergence
%
%  Plot the convergence results of WMRNSD, KWMRNSD, and unconstrained LSQR
%  for the test problems used in IR_ConstrainedExamples script

axes('FontSize',FS)
plot( IterInfo_lsqr.Enrm,'-bd','LineWidth',LW,'MarkerSize',MS ), hold on
xlim([0 MaxIter+5])
plot( IterInfo_wmrnsd.Enrm,'-m^','LineWidth',LW,'MarkerSize',MS)
plot( IterInfo_kwmrnsd.Enrm,'-cV','LineWidth',LW,'MarkerSize',MS)
legend('LSQR', 'WMRNSD', 'KWMRNSD','Location','Best');
delete(findobj(gca,'Type','line','Marker','d'));
delete(findobj(gca,'Type','line','Marker','^'));
delete(findobj(gca,'Type','line','Marker','V'));
plot( IterInfo_lsqr.Enrm,'b','LineWidth',LW )
plot( IterInfo_wmrnsd.Enrm,'m','LineWidth',LW)
plot( IterInfo_kwmrnsd.Enrm,'c','LineWidth',LW)
plot(1:11:MaxIter,IterInfo_lsqr.Enrm(1:11:MaxIter),'bd','MarkerSize',MS, 'LineWidth',LW)
plot(1:11:MaxIter,IterInfo_wmrnsd.Enrm(1:11:MaxIter),'m^','MarkerSize',MS,'LineWidth',LW)
plot(1:11:MaxIter,IterInfo_kwmrnsd.Enrm(1:11:MaxIter),'cV','MarkerSize',MS,'LineWidth',LW)
xlabel('Iterations','fontsize',FS)
ylabel('Relative error','fontsize',FS)