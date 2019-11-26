function [] = population_plot(t, sol, plottype, save, name)
% Plot populations of T cells, Natural Killer cells, and tumor cells over
% time and save if desired

if plottype == 0 % separate plots
    subplot(3,1,1);
    plot(t, sol(:, 1));
    xlabel('Days');
    ylabel('Cell population');
    title('Tumor cells');
    
    subplot(3,1,2);
    plot(t, sol(:, 2));
    xlabel('Days');
    ylabel('Cell population');
    title('Natural Killer cells');
    
    subplot(3,1,3);
    plot(t, sol(:, 3));
    xlabel('Days');
    ylabel('Cell population');
    title('CD8+ T cells');
    
elseif plottype == 1 % plot all together
    plot(t, sol);
    set(gca,'Yscale','log')
    xlabel('Days');
    ylabel('Cell population');
    legend(["T", "N", "L"]);
    title("One nice lookin' plot");
end
if save
   savefig(name);
end