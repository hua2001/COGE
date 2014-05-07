function [] = EnergyPlot( )
%=========================================================================

load('energy.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure1 = figure('NumberTitle','off','Name','TotalEnergy');
% Create axes
axes1 = axes('Parent',figure1,'FontSize',18);
xlim(axes1,[0 1100]);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
% Create multiple lines using matrix input to plot
plot1 = plot(r_energy_tensile(:,1),'Parent',axes1,'LineWidth',3);
plot2 = plot(r_energy_shear(:,1),'Parent',axes1,'Linewidth',3);
plot3 = plot(r_energy_compression(:,1),'Parent',axes1,'Linewidth',3);
%%%%%
set(plot1,'Color',[0 0 1],'DisplayName','Uniaxial Tensile Test');
set(plot2,'LineStyle','--','Color',[1 0 0],...
    'DisplayName','Shear Test');
set(plot3,'LineStyle','-.','Color',[0 1 1],...
    'DisplayName','Triaxial Compression Test');
% Create xlabel
xlabel('Substeps','FontSize',20);
% Create ylabel
ylabel('Energy','FontSize',20);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','NorthWest','FontSize',15);
print -depsc 'TotalEnergy.eps'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure1 = figure('NumberTitle','off','Name','ElasticStrainEnergy');
% Create axes
axes1 = axes('Parent',figure1,'FontSize',18);
xlim(axes1,[0 1100]);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
% Create multiple lines using matrix input to plot
plot1 = plot(r_energy_tensile(:,7),'Parent',axes1,'LineWidth',3);
plot2 = plot(r_energy_shear(:,7),'Parent',axes1,'Linewidth',3);
plot3 = plot(r_energy_compression(:,7),'Parent',axes1,'Linewidth',3);
%
set(plot1,'Color',[0 0 1],'DisplayName','Uniaxial Tensile Test');
set(plot2,'LineStyle','--','Color',[1 0 0],...
    'DisplayName','Shear Test');
set(plot3,'LineStyle','-.','Color',[0 1 1],...
    'DisplayName','Triaxial Compression Test');
% Create xlabel
xlabel('Substeps','FontSize',20);
% Create ylabel
ylabel('Energy','FontSize',20);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','NorthWest','FontSize',15);
print -depsc 'ElasticStrainEnergy.eps'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure1 = figure('NumberTitle','off','Name','IrreveasibleStrainEnergy');
% Create axes
axes1 = axes('Parent',figure1,'FontSize',18);
xlim(axes1,[0 1100]);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
% Create multiple lines using matrix input to plot
plot1 = plot(r_energy_tensile(:,4),'Parent',axes1,'LineWidth',3);
plot2 = plot(r_energy_shear(:,4),'Parent',axes1,'Linewidth',3);
plot3 = plot(r_energy_compression(:,4),'Parent',axes1,'Linewidth',3);
%%%%%
set(plot1,'Color',[0 0 1],'DisplayName','Uniaxial Tensile Test');
set(plot2,'LineStyle','--','Color',[1 0 0],...
    'DisplayName','Shear Test');
set(plot3,'LineStyle','-.','Color',[0 1 1],...
    'DisplayName','Triaxial Compression Test');
% Create xlabel
xlabel('Substeps','FontSize',20);
% Create ylabel
ylabel('Energy','FontSize',20);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','NorthWest','FontSize',15);
print -depsc 'IrreveasibleStrainEnergy.eps'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure1 = figure('NumberTitle','off','Name','CrackDebondingEnergy');
% Create axes
axes1 = axes('Parent',figure1,'FontSize',18);
xlim(axes1,[0 1100]);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
% Create multiple lines using matrix input to plot
plot1 = plot(r_energy_tensile(:,5),'Parent',axes1,'LineWidth',3);
plot2 = plot(r_energy_shear(:,5),'Parent',axes1,'Linewidth',3);
plot3 = plot(r_energy_compression(:,5),'Parent',axes1,'Linewidth',3);
%%%%%
set(plot1,'Color',[0 0 1],'DisplayName','Uniaxial Tensile Test');
set(plot2,'LineStyle','--','Color',[1 0 0],...
    'DisplayName','Shear Test');
set(plot3,'LineStyle','-.','Color',[0 1 1],...
    'DisplayName','Triaxial Compression Test');
% Create xlabel
xlabel('Substeps','FontSize',20);
% Create ylabel
ylabel('Energy','FontSize',20);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','NorthWest','FontSize',15);
print -depsc 'CrackDebondingEnergy.eps'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure1 = figure('NumberTitle','off','Name','ElasticStrainEnergyPercentage');
% Create axes
axes1 = axes('Parent',figure1,'FontSize',18);
xlim(axes1,[0 1100]);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
% Create multiple lines using matrix input to plot
plot1 = plot(r_energy_tensile(:,7)./r_energy_tensile(:,1),'Parent',axes1,'LineWidth',3);
plot2 = plot(r_energy_shear(:,7)./r_energy_shear(:,1),'Parent',axes1,'Linewidth',3);
plot3 = plot(r_energy_compression(:,7)./r_energy_compression(:,1),'Parent',axes1,'Linewidth',3);
%
set(plot1,'Color',[0 0 1],'DisplayName','Uniaxial Tensile Test');
set(plot2,'LineStyle','--','Color',[1 0 0],...
    'DisplayName','Shear Test');
set(plot3,'LineStyle','-.','Color',[0 1 1],...
    'DisplayName','Triaxial Compression Test');
% Create xlabel
xlabel('Substeps','FontSize',20);
% Create ylabel
ylabel('Energy','FontSize',20);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','Best','FontSize',15);
print -depsc 'ElasticStrainEnergyPercentage.eps'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure1 = figure('NumberTitle','off','Name','IrreveasibleStrainEnergyPercentage');
% Create axes
axes1 = axes('Parent',figure1,'FontSize',18);
xlim(axes1,[0 1100]);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
% Create multiple lines using matrix input to plot
plot1 = plot(r_energy_tensile(:,4)./r_energy_tensile(:,1),'Parent',axes1,'LineWidth',3);
plot2 = plot(r_energy_shear(:,4)./r_energy_shear(:,1),'Parent',axes1,'Linewidth',3);
plot3 = plot(r_energy_compression(:,4)./r_energy_compression(:,1),'Parent',axes1,'Linewidth',3);
%%%%%
set(plot1,'Color',[0 0 1],'DisplayName','Uniaxial Tensile Test');
set(plot2,'LineStyle','--','Color',[1 0 0],...
    'DisplayName','Shear Test');
set(plot3,'LineStyle','-.','Color',[0 1 1],...
    'DisplayName','Triaxial Compression Test');
% Create xlabel
xlabel('Substeps','FontSize',20);
% Create ylabel
ylabel('Energy','FontSize',20);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','NorthWest','FontSize',15);
print -depsc 'IrreveasibleStrainEnergyPercentage.eps'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure1 = figure('NumberTitle','off','Name','CrackDebondingEnergyPercentage');
% Create axes
axes1 = axes('Parent',figure1,'FontSize',18);
xlim(axes1,[0 1100]);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
% Create multiple lines using matrix input to plot
plot1 = plot(r_energy_tensile(:,5)./r_energy_tensile(:,1),'Parent',axes1,'LineWidth',3);
plot2 = plot(r_energy_shear(:,5)./r_energy_shear(:,1),'Parent',axes1,'Linewidth',3);
plot3 = plot(r_energy_compression(:,5)./r_energy_compression(:,1),'Parent',axes1,'Linewidth',3);
%%%%%
set(plot1,'Color',[0 0 1],'DisplayName','Uniaxial Tensile Test');
set(plot2,'LineStyle','--','Color',[1 0 0],...
    'DisplayName','Shear Test');
set(plot3,'LineStyle','-.','Color',[0 1 1],...
    'DisplayName','Triaxial Compression Test');
% Create xlabel
xlabel('Substeps','FontSize',20);
% Create ylabel
ylabel('Energy','FontSize',20);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','NorthWest','FontSize',15);
print -depsc 'CrackDebondingEnergyPercentage.eps'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure1 = figure('NumberTitle','off','Name','EnergyPercentage');
% Create axes
axes1 = axes('Parent',figure1,'FontSize',18);
xlim(axes1,[0 1100]);
ylim(axes1,[0 100]);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
% Create multiple lines using matrix input to plot
end_r = length(r_energy);
plot1 = plot(100.*r_energy(2:end_r,1)./r_energy(2:end_r,1),'Parent',axes1,'LineWidth',3);
plot2 = plot(100.*r_energy(2:end_r,7)./r_energy(2:end_r,1),'Parent',axes1,'Linewidth',3);
plot3 = plot(100.*r_energy(2:end_r,4)./r_energy(2:end_r,1),'Parent',axes1,'Linewidth',5);
plot4 = plot(100.*r_energy(2:end_r,5)./r_energy(2:end_r,1),'Parent',axes1,'Linewidth',3);

set(plot1,'Color',[0 0 1],'DisplayName','External work');
set(plot2,'LineStyle','--','Color',[1 0 0],...
    'DisplayName','Elastic strain energy');
set(plot3,'LineStyle',':','Color',[1 0.5 0],...
    'DisplayName','Irreveasible strain energy');
set(plot4,'LineStyle','-.','Color',[0 1 1],...
    'DisplayName','Crack debonding');
% Create xlabel
xlabel('Substeps','FontSize',20);
% Create ylabel
ylabel('E_i/E_{Total}, %','FontSize',20);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','NorthWest','FontSize',15);
print -depsc 'energypercentage.eps'
end


end

