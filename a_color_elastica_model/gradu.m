function [ux,uy]=gradu(u)

uex=expandf3(u);
ux=(uex(3:end,2:end-1,:)-uex(2:end-1,2:end-1,:));
uy=(uex(2:end-1,3:end,:)-uex(2:end-1,2:end-1,:));