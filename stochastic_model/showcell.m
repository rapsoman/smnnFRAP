function showcell(YG1,XG1,YG2,ZG2,YF1,XF1,YF2,ZF2,q,x,isbleached,showbleach)

subplot(2,1,2);
plot(YG1,-XG1,'b.');
axis([-5,5,-5,5]);
hold on
plot(YG1,-XG1,'m.');
if showbleach
    plot(YF1,-XF1,'k.');
end

subplot(2,1,1);
plot(YG2,ZG2,'b.');
axis([-5,5,-5,5]);
hold on
plot(YG2,ZG2,'m.');
if showbleach
    plot(YF2,ZF2,'k.');
end

N=size(q,1);
n=3;

bleached_bound=find(isbleached & q);
notbleached_bound=find((~isbleached) & q);
bleached_notbound=find(isbleached & (~q));
notbleached_notbound=find((~isbleached) & (~q));

x1=(0:n:n*(N-1))+1;
x2=(0:n:n*(N-1))+2;
x3=(0:n:n*(N-1))+3;

subplot(2,1,2);


if ~isempty(notbleached_notbound), plot(x(x2(notbleached_notbound)),-x(x1(notbleached_notbound)),'.g'); end
if ~isempty(notbleached_bound), plot(x(x2(notbleached_bound)),-x(x1(notbleached_bound)),'xg'); end
if ~isempty(bleached_notbound), plot(x(x2(bleached_notbound)),-x(x1(bleached_notbound)),'.k'); end
if ~isempty(bleached_bound), plot(x(x2(bleached_bound)),-x(x1(bleached_bound)),'xk'); end


subplot(2,1,1);


if ~isempty(notbleached_notbound), plot(x(x2(notbleached_notbound)),x(x3(notbleached_notbound)),'.g'); end
if ~isempty(notbleached_bound), plot(x(x2(notbleached_bound)),x(x3(notbleached_bound)),'xg'); end
if ~isempty(bleached_notbound), plot(x(x2(bleached_notbound)),x(x3(bleached_notbound)),'.k'); end
if ~isempty(bleached_bound), plot(x(x2(bleached_bound)),x(x3(bleached_bound)),'xk'); end



drawnow;


subplot(2,1,1);
hold off

subplot(2,1,2);
hold off