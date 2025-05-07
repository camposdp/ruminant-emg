function [] = getSigstar(b,es,p,pcrit)

 for i = 1:length(p)
     if p(i) < pcrit
         y = es(i,:)*100;
         ctr1 = bsxfun(@plus, b(1).XData, [b(1).XOffset]');
         ctr2 = bsxfun(@plus, b(2).XData, [b(2).XOffset]');
         ctr = linspace(ctr1(i),ctr2(i),2);
         plot(ctr,...
             [1 1]*max(y)+1, '-k', 'LineWidth',1)
         center = mean([ctr2(i) ctr1(i)]);
         
         v1 = linspace(max(y)+0.5,max(y)+1,2);
         plot([1 1]*ctr1(i),v1,'-k', 'LineWidth',1)
         
         v2 = linspace(max(y)+0.5,max(y)+1,2);
         plot([1 1]*ctr2(i),v2,'-k', 'LineWidth',1)
         
         
         plot(center,max(y)+2, '*k')
     end
 end
 
end