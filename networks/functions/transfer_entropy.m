function [T,H,S] = transfer_entropy(Ys,Xt,Yt,dx,dy)

[N,D] = size(Xt); assert(N == size(Yt,1)); assert(D==1);
[N,D] = size(Yt); assert(N == size(Ys,1)); assert(D==1);
[N,D] = size(Ys); assert(N == size(Xt,1)); assert(D==1);

assert(isempty(find(isnan(Ys),1)));
assert(isempty(find(isnan(Xt),1)));
assert(isempty(find(isnan(Yt),1)));

Bmin = min(Xt)-1e-6; 
Bmax = max(Xt)+1e-6; 
Bx = Bmin:dx:Bmax;

Bmin = min(min(Yt),min(Ys))-1e-6; 
Bmax = max(max(Yt),max(Ys))+1e-6; 
By = Bmin:dy:Bmax;

[Psxt,Psx,Pst,Pxt,Ps,Px,Pt] = hist3(Ys,Xt,Yt,By,Bx,By);

Psxt = Psxt(:);
Psx  = Psx(:);
Pst  = Pst(:);
Pxt  = Pxt(:);
Ps   = Ps(:);
Px   = Px(:);
Pt   = Pt(:);

if abs(sum(Ps)-1)>1/N^2; error(' ');end;%('Ps does not sum to 1'); end;
if abs(sum(Px)-1)>1/N^2; error(' ');end;%('Px does not sum to 1'); end;
if abs(sum(Pt)-1)>1/N^2; error(' ');end;%('Pt does not sum to 1'); end;

if abs(sum(Psx)-1)>1/N^2; error(' ');end;%('Psx does not sum to 1'); end;
if abs(sum(Pst)-1)>1/N^2; error(' ');end;%('Pst does not sum to 1'); end;
if abs(sum(Pxt)-1)>1/N^2; error(' ');end;%('Pxt does not sum to 1'); end;

if abs(sum(Psxt)-1)>1/N^2; error(' ');end;%('Psxt does not sum to 1'); end;

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






