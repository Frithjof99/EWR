domain = 'S';
geom = [1/2 1/2 1 1];
NperDim = 17;
star = 5;

order1 = 'central';

fun_a = 1;
fun_b = [];
fun_c = [];

fun_f = [];

fun_uex = [];

fun_dir = @(x,y)0;

fun_nlf = @(u)lambda*exp(u);
fun_nlfp = @(u)lambda*exp(u);

gopt = fdplot();
gopt.zlbl = 'u_h(x,y)';
gopt.fixaxis = [0 1 0 1 0 1.3];