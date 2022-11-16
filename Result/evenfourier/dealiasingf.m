function A = dealiasingf(A)

A(floor(end/3)+1:end-floor(end/3)+1,:)=0;
A(:,floor(end/3)+1:end-floor(end/3)+1)=0;
end

