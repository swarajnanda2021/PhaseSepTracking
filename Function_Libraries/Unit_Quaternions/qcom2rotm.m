function R = qcom2rotm( q )
%qver2rotm convert unit quaternion, i.e. versor, to a rotationmatrix in 3D
%space, suitable for symbolic usage. Koen Muller

% unit normalization sq
nsq=dot(q,q);

% get components
q0=q(1);
qi=q(2);
qj=q(3);
qk=q(4); % formalism st q_ver=(q0,q_vec)^T, note that q_ver contains components, and does not operate as a vector.

% define complex part
qc=[qi
    qj
    qk];

% define cross operator
Q=[0 -qk qj
    qk 0 -qi
    -qj qi 0];

% rotation matrix as define in formalism
R=(q0^2-qc'*qc)*eye(3)+2*(qc*qc')+2*q0*Q;

% normalize such that input not need to be unitary?
R=R/nsq; % @ norm, q is redefined

end