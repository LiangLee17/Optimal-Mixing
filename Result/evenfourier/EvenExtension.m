function A = EvenExtension(A,ec)
%output the even extension of A.
%or reverse the even extension. 


if ec=='e'
  A = [A; flipud(A(2:end-1,:))];
A = [A fliplr(A(:,2:end-1))];  
return
end
if ec=='c'
    A=A(1:end/2+1,1:end/2+1);
    return
end
end

