clear StatenewViterfind StatenewViter Delta Shi
        for i=1:N
            Delta(1,i)=Pi(i)*b(i,Ob(1));
            Shi(1,i)=0;
        end
        for t=2:T
            for j=1:N
                for i=1:N
                    tofindmaxshi(i)=Delta(t-1,i)*a(i,j);
                    tofindmaxdel(i)=Delta(t-1,i)*a(i,j);
                end
                Delta(t,j)=max(tofindmaxdel)*b(j,Ob(t));
               [temp, Shi(t,j)]=max(tofindmaxshi);
%                 Delta(t,j)=Delta(t,j)*b(j,Ob(t));
                clear tofindmax;
            end
        end
        for i=1:N
            tofindmax(i)=Delta(T,i);
        end
        [Pnew,StatenewViter(T)]=max(tofindmax);
        
        for t=T-1:-1:1
            StatenewViter(t)=Shi(t+1,StatenewViter(t+1));
        end
        