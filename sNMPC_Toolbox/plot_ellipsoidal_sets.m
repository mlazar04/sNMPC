%% Plot the obtained ellipsoidal terminal sets
function plot_ellipsoidal_sets(sys, p, E, VOL1, XUset, Xset_scaled)

try
    testET=ellipsoid(eye(2));
    %testET.getShapeMat; % remove this
    ET=1;
catch
    ET=0;
end

n=sys.n; x_low=sys.x_low; x_high=sys.x_high;
M=p.M;

[maxvol,maxvol_ind]=max(cell2mat(VOL1));
if n<3
    %figure
    hold on
    for i=2:M-1
        if ET==1
            plot(E{1,i},'b');
        else
            plot(E{1,i},'alpha',0,'edgecolor','blue');
        end
    end
    
    if ET==1
        plot(E{1,1},'b');
    else
        plot(E{1,1},'alpha',0,'edgecolor','blue');
        plot(E{1,end},'alpha',0,'edgecolor','blue');
    end
    
    % Plot state constraints and set of admissable sates corresponding to
    % the largest set
    plot(Xset_scaled,'alpha',0.1);
    pause(0.01);
else
    % Only possible with the Ellipsoidal Toolbox (for higher dim. systems)
    for i=1:2:n
        %figure()
        hold on
        for o=1:M
            hold on
            if rem(n,2)==1 && i==n
                projE=projection(E{1,o}, [[zeros(i-2,1);1;zeros(n-i+1,1)] [zeros(i-1,1);1;zeros(n-1-i,1)]]);
            else
                projE=projection(E{1,o}, [[zeros(i-1,1);1;zeros(n-i,1)] [zeros(i,1);1;zeros(n-1-i,1)]]);
            end
            if o==1
                plot(projE,'b'); 
            elseif o==M
                plot(projE,'b'); 
            else
                plot(projE,'b'); 
            end
        end
        if rem(n,2)==0
            xlim([x_low(i) x_high(i)])
            ylim([x_low(i+1) x_high(i+1)])
        end
    end
    pause(0.3);
end
