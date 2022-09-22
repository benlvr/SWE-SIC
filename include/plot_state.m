function plot_state(xc,xe,t,t_idx,eta_old,u_old,Hu,eta_axis,u_axis,q_axis,stride)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if mod(t_idx,stride) == 0
    drawnow
    subplot(1,3,1)
    plot(xc,eta_old(2:end-1))
    axis(eta_axis)
    xlabel('x')
    ylabel('\eta')
    subplot(1,3,2)
    plot(xe,u_old(2:end-1))
    axis(u_axis)
    ylabel('u')
    subplot(1,3,3)
    plot(xe,Hu.*u_old(2:end-1))
    axis(q_axis)
    ylabel('q')
    title(t)
end

end

