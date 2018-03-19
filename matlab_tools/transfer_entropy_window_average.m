function [T,H,S] = transfer_entropy_window_average(Ys,Xt,Yt,Bx,By)

[N,D] = size(Xt);
assert(N==size(Yt,1));assert(D==size(Yt,2));
assert(N==size(Ys,1));assert(D==size(Ys,2));

Mx = length(Bx);
My = length(By);

Psxt = zeros(My,Mx,My);

for d = 1:D

 for s = 2:My
  Is = find(Ys(:,d) < By(s));
  Iss = find(Ys(Is,d) >= By(s-1));
  Is = Is(Iss);
 
  if ~isempty(Is)
   for x = 2:Mx
    Ix = find(Xt(Is) < Bx(x));
    Ixx = find(Xt(Is(Ix)) >= Bx(x-1));
    Isx = Is(Ix(Ixx));
  
    if ~isempty(Isx)
     for t = 2:My
      It = find(Yt(Isx) < By(t));
      Itt = find(Yt(Isx(It)) >= By(t-1));
      Isxt = Isx(It(Itt));
   
      if ~isempty(Isxt)
       Psxt(s,x,t) = Psxt(s,x,t) + length(Isxt);
      end % if Isxt
     end % t
 
    end % if Isx
   end % x
  
  end % if Is
 end % s

end %d
 
N = N*D;
Psxt = Psxt/N;

Psx = squeeze(sum(Psxt,3));
Pst = squeeze(sum(Psxt,2));
Pxt = squeeze(sum(Psxt,1));
Ps  = squeeze(sum(Psx,2));
Px  = squeeze(sum(Psx,1));
Pt  = squeeze(sum(Pst,1));

Ps   = Ps(:);
Px   = Px(:);
Pt   = Pt(:);
Pxt  = Pxt(:);
Pst  = Pst(:);
Psx  = Psx(:);
Psxt = Psxt(:);

if abs(sum(Ps)-1)>1/N^2; keyboard;end;%('Ps does not sum to 1'); end;
if abs(sum(Px)-1)>1/N^2; keyboard;end;%('Px does not sum to 1'); end;
if abs(sum(Pt)-1)>1/N^2; keyboard;end;%('Pt does not sum to 1'); end;

if abs(sum(Psx)-1)>1/N^2; keyboard;end;%('Psx does not sum to 1'); end;
if abs(sum(Pst)-1)>1/N^2; keyboard;end;%('Pst does not sum to 1'); end;
if abs(sum(Pxt)-1)>1/N^2; keyboard;end;%('Pxt does not sum to 1'); end;

if abs(sum(Psxt)-1)>1/N^2; keyboard;end;%('Psxt does not sum to 1'); end;

Hs = -Ps(Ps>0)'*log(Ps(Ps>0));
Hx = -Px(Px>0)'*log(Px(Px>0));
Ht = -Pt(Pt>0)'*log(Pt(Pt>0));

Hsx = -Psx(Psx>0)'*log(Psx(Psx>0));
Hst = -Pst(Pst>0)'*log(Pst(Pst>0));
Hxt = -Pxt(Pxt>0)'*log(Pxt(Pxt>0));

Hsxt = -Psxt(Psxt>0)'*log(Psxt(Psxt>0));

Isx = Hs+Hx-Hsx;
Ist = Hs+Ht-Hst;
Ixt = Hx+Ht-Hxt;
Isxt = -(Hs+Hx+Ht-Hsxt-Isx-Ist-Ixt);

T = (Isx-Isxt);
H = Hs;
S = Isxt;






