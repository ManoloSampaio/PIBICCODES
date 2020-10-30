function boleano=isstable(A)
    V =eig(A);
    for k =1:length(V)
        if(V(k)>0)
            boleano = 0; %The system is not stable.
            return;
        end
        if(V(k)==0)
            boleano = 2;
            return;%The system is marginaly stable.
        end
    end
    boleano = 3;
end