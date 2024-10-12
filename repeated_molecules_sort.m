% RECURSIVE FUNCTION SORTS MOLECULES ------------------------- 

function [final_two,g_sort]=repeated_molecules_sort(i,k,b,final_two,g_sort)
%g_sort(i,1:2)=final_two(i,1:2);
g_sort(i,k)=final_two(b,2);
a1=final_two(b,2);
final_two(b,1:2)=0;
[x2 y2]=find(final_two==a1);
    if (size(x2,1)~=0)
      k=k+1;
      c=x2(1,1);
      [final_two,g_sort]=repeated_molecules_sort(i,k,c,final_two,g_sort);
    end
end

%---------------------------------------------------------------

