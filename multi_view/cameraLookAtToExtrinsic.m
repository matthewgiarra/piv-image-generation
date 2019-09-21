

function [M, R, t] = cameraLookAtToExtrinsic(EYE, CENTER, UP)

    % OpenGL - like "eye, center, up" camera matrix
    % Details: http://ksimek.github.io/2012/08/22/extrinsic/
    
    L = CENTER-EYE;
    Ln = L ./ norm(L,2);
    s = cross(Ln, UP);
    sn = s ./ norm(s,2);
    uprime = cross(sn,Ln);
    R = [sn; uprime; -Ln];
    t = -R*EYE';
    M = [R, t; 0,0,0,1];

end
