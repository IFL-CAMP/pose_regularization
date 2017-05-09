function M = genPose(alpha,beta,gamma,x,y,z)

M = ...
[                                    cos(beta)*cos(gamma),                                   -cos(beta)*sin(gamma),            -sin(beta), x;...
  cos(alpha)*sin(gamma) - cos(gamma)*sin(alpha)*sin(beta), cos(alpha)*cos(gamma) + sin(alpha)*sin(beta)*sin(gamma), -cos(beta)*sin(alpha), y;...
  sin(alpha)*sin(gamma) + cos(alpha)*cos(gamma)*sin(beta), cos(gamma)*sin(alpha) - cos(alpha)*sin(beta)*sin(gamma),  cos(alpha)*cos(beta), z;...
                                                        0,                                                       0,                     0, 1];
                                                    
end