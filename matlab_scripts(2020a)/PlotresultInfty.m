% This script plots all the figures for the infty. 
close all
Vec = [1, 6, 11, 26, 56, length(A)];
for k=1:6
    iteration = Vec(k);
    f=figure;
    f.Position = [680 753 334 345];
    f.Color = 'W';
    plot3(x0,y0,z0,'-b','LineWidth', 4)
    hold all
    plot3(ones(size(x0))*4,y0,z0,'-k', 'LineWidth', 2.5)
    plot3(x0,ones(size(y0))*4.001,z0,'-k', 'LineWidth', 2.5)
    plot3(x0,y0,ones(size(z0))*(-2),'-k', 'LineWidth', 2.5)
    X_stars = A{iteration};
    X = splinepoints(A{iteration},11);
    plot3(X(1,:),X(2,:),X(3,:),'-','LineWidth', 1,'Color','r')
    hold all
    plot3(X_stars(1,:),X_stars(2,:),X_stars(3,:), '*','MarkerSize',3.5, 'LineWidth', 1,'Color','r')
    plot3(ones(size(X(1,:)))*4,X(2,:),X(3,:),'-r', 'LineWidth', 2.5)
    plot3(X(1,:),ones(size(X(2,:)))*4.001,X(3,:),'-r', 'LineWidth', 2.5)
    plot3(X(1,:),X(2,:),ones(size(X(3,:)))*(-2),'-r', 'LineWidth', 2.5)
    hold off
    axis equal
    axis([-4.,4.,-4.001,4.001,-2,6])
    set(gca,'fontsize',14)
    grid on
    set(gca,'GridAlpha', 0.2);
    set(gca,'LineWidth',2.,'TickLength',[0.025 0.04]);
    hold off
    if k==1
        title(strcat("$\ell=$",num2str(iteration-1)," (initial guess)"),'Interpreter','Latex','FontSize',20)
    elseif k==6
        title(strcat("$\ell=$",num2str(iteration-1)," (final result)"),'Interpreter','Latex','FontSize',20)
    else
        title(strcat("$\ell=$",num2str(iteration-1)),'Interpreter','Latex','FontSize',20)
    end
    if save_results == 1
        filename=strcat(num2str(k),name);
        print(gcf,'-djpeg',filename);
        print(gcf,'-depsc',filename);
        savefig(filename);
    end
end