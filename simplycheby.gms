
$title example gams code of Chebyshev and ordinary polynomial regression.

$include cheby.gms

option decimals=8, optcr=0.0,reslim=9E9,threads=-1;

*1=ordinary polynomial; 2=chebyshev polynomial
$set type 2

sets i_x /0*10/, i_b /0*5/;

*scalar min_x /-1/, max_x /3/, delta /0/;
scalar min_x /-5/, max_x /5/, delta /0/;

parameter x(i_x),wgt_x(i_x);

if(%type% = 1,
*ORDINARY POLYNOMIAL
        delta = (max_x - min_x)/(card(i_x)-1);
        x(i_x) = min_x + delta*i_x.val;
        wgt_x(i_x)=1;
)

if(%type% = 2,
*CHEBYSHEV POLYNOMIAL
        x(i_x) = shifstategri(chebycolnod(i_x),min_x,max_x);
        wgt_x(i_x) = (1-sqr(chebycolnod(i_x)))**(-1/2);
);

display x, wgt_x;
*$exit;

$macro f(x) (1/3)*power(x,3) + 2*power(x,2) + x - 10
*$macro f(x) 1/(1+power(x,2))

parameter y(i_x);

y(i_x) = f(x(i_x));

*display y;
*$exit

variables coef(i_b), error(i_x),obj,yapprox(i_x);

equations eqn_error(i_x),eqn_obj,eqnyapprox(i_x);

eqnyapprox(i_x)..       yapprox(i_x) =e= sum(i_b, coef(i_b)*power(x(i_x),(ord(i_b)-1)))$(%type% = 1)
                                                                + sum(i_b, coef(i_b)*(chebyrecgen(i_b,(normstategri(x(i_x),min_x,max_x)))))$(%type% = 2);
eqn_error(i_x)..        error(i_x) =e= y(i_x) - yapprox(i_x);
eqn_obj..       obj =e= sum(i_x,wgt_x(i_x)*error(i_x)*error(i_x));

model cheby /eqn_obj, eqn_error,eqnyapprox/;

option nlp = ipopt;
solve cheby using nlp minimizing obj;

*display x, y,yapprox.l,error.l,coef.l;

sets j /y,ya/;
parameter test(i_x,*);
test(i_x,'y(x)') = y(i_x);
test(i_x,'y(x) approx') = yapprox.l(i_x);

display test;
*display y,yapprox.l;
*$exit;

*$libinclude gnuplotxyz test
*$libinclude gnuplotxyz y x-axis y-axis;
*$libinclude gnuplotxyz yapprox.l x-axis y-axis;
$exit;


