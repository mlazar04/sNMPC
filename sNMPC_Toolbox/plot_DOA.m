%% Plot domain of attraction representation of the closed-loop system return number of feasible and infeasible gridpoints
function [feasible_points, infeasible_points] = plot_DOA(s,sys,p,P,alpha,gammascale,x_axis,y_axis,E2, VOL2, XUset, Xset_scaled)

feasible_points = 0;
infeasible_points = 0;

figure(); hold on;
xlabel('x1') 
ylabel('x2') 
plot_ellipsoidal_sets(sys, p, E2, VOL2, XUset, Xset_scaled);
pause(0.5);

for l=1:1:length(x_axis)
    for j=1:1:length(y_axis)
        
        % Check if problem is feasible for any set of terminal ingredients
            x0=[x_axis(l);y_axis(j)];
            [feasible,init_index]=find_init_set(s,p,P,alpha,gammascale,x0); % See if an initial set can be found
         if feasible
            plot(x_axis(l),y_axis(j),'o','MarkerFaceColor','g')
            feasible_points = feasible_points+1;
         else % If problem not feasible
            plot(x_axis(l),y_axis(j),'o','MarkerFaceColor','r')
            infeasible_points = infeasible_points+1;
        end
    end
end
end