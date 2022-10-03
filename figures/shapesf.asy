settings.outformat="pdf";
texpreamble("\usepackage[charter]{mathdesign}");
//defaultpen(fontsize(1pt));

unitsize(4cm);
dotfactor=6;

real width_box=2.0;
real height_box=3.0;
real arrow_size=10;
real size_circle=0.09;
real lw=2;
int n=20;


//---------------Definitions of functions start-------------------------------------------------------------------------------
pair find_midpoint(pair a, pair b){
	pair c = (a.x+b.x, a.y+b.y);
	return (c.x*0.5 , c.y*0.5);	
}

void draw_rect(pair x, real w=width_box, real h=height_box, real lw=1.0){
	draw(x--(x.x+w, x.y)--(x.x+w, x.y+h)--(x.x,x.y+h)--cycle, black+linewidth(lw));
}

void draw_tri(pair lb=(0,0), real size=1.0){
	pair b=lb;
	pair c=(lb.x+size, lb.y);
	pair a=extension(b,b+dir(60),c,c+dir(120));
	//pair d=extension(b,b+dir(30),a,a+dir(270));
	draw(a--b--c--cycle, black+linewidth(lw));

	pair tri_vert[]={a,b,c};

	//EDIT VERTEX LABELS OF TRIANGLE
	label("$\textbf{GRE2}$", a, 4N, black+fontsize(n));
	label("HEL2", b, SW, black+fontsize(n));
	label("MGM1", c, SE, black+fontsize(n));

	label("$\textbf{0.13}$", a, 0.5S+1.5N, brown+fontsize(n));
	label("-1.06", b, 4S+5W, brown+fontsize(n));
	label("0.05", c, 4S+4E, brown+fontsize(n));

	//EDIT EDGE LABELS OF TRIANGLE
	label("-.0283", find_midpoint(a,b), NW, blue+fontsize(n));
	label("0.2246", find_midpoint(b,c), S, blue+fontsize(n));
	label("-0.1528", find_midpoint(c,a), NE, blue+fontsize(n));


	for(int i=0;i<3;++i){
		//filldraw(shift(tri_vert[i])*scale(size_circle)*unitcircle, white, black+linewidth(width_line));
	}
}

//draw_tri();


void draw_sq(pair ld=(0,0), real scale=1.0 ){
	real size=1.414*scale;
	pair b=ld;
	pair c=(ld.x+size, ld.y);
	pair d=extension(b,b+dir(45),c,c+dir(135));
	pair e=extension(b,b+dir(315),c,c+dir(225));
	draw(b--e--c--d--cycle, black+linewidth(lw));

	pair sq_vert[]={b,c,d,e};
	
	//EDIT VERTEX LABELS OF SQUARE
	label("HEL2", b, W, black+fontsize(n));
	label("DPH2", c, E, black+fontsize(n));
	label("$\textbf{GRE2}$", d, 4N, black+fontsize(n));
	label("ATM1", e, S, black+fontsize(n));

	label("-0.27", b, 2.5S+4W, brown+fontsize(n));
	label("-0.09", c, 2.5S+3.5E, brown+fontsize(n));
	label("$\textbf{-0.47}$", d, 0.5S+1.5N,brown+fontsize(n));
	label("0.45", e, 4S, brown+fontsize(n));

	//EDIT EDGE LABELS OF SQUARE
	label("", find_midpoint(b,e), SW, blue+fontsize(n));
	label("NA", find_midpoint(c,e), SE, blue+fontsize(n));
	label("-0.2579", find_midpoint(d,b), NW, blue+fontsize(n));
	label("-0.03135", find_midpoint(c,d), NE, blue+fontsize(n));

	real side=size/1.414;
	pair ln=(b.x-0.5*side,b.y);
	pair lf=(b.x-0.5*side*2,b.y);
	pair rn=(c.x+0.5*side,c.y);
	pair rf=(c.x+0.5*side*2,c.y);

	//draw(e--ln--d, black+dotted);
	//draw(e--lf--d, black+dotted);
	//draw(e--rn--d, black+dotted);
	//draw(e--rf--d, black+dotted);

	for(int i=0;i<4;++i){
		//filldraw(shift(sq_vert[i])*scale(size_circle)*unitcircle, white, black+linewidth(width_line));
	}
	 
}

void draw_pent(pair lb=(0,0), real size=1.0){
	pair a=lb;
	pair b=(lb.x+size, lb.y);
	pair c=extension(a,a+dir(36),b,b+dir(72));
	pair d=extension(a,a+dir(72),b,b+dir(108));
	pair e=extension(a,a+dir(108),b,b+dir(144));
	draw(a--b--c--d--e--cycle, black+linewidth(lw));

	pair pent_vert[]={a,b,c,d,e};
	
	//EDIT VERTEX LABELS OF PENTAGON HERE 
	label("YBL071W", a, W+S, black+fontsize(n));
	label("YGR195W", b, E+S, black+fontsize(n));
	label("YGL078C", c, E, black+fontsize(n));
	label("$\textbf{YLR120C}$", d, 4N, black+fontsize(n));
	label("YFR033C", e, W, black+fontsize(n));

	label("-0.18", a, 4S+4W, brown+fontsize(n));
	label("0.12", b, 4S+4E, brown+fontsize(n));
	label("0.07", c, 2.5S+4E, brown+fontsize(n));
	label("$\textbf{0.15}$", d, N, brown+fontsize(n));
	label("0.31", e, 2.5S+4W, brown+fontsize(n));

	//EDIT EDGE LABELS OF PENTAGON HERE
	label("0.0052", find_midpoint(a,b), S, blue+fontsize(n));
	label("NA", find_midpoint(b,c), SE, blue+fontsize(n));
	label("NA", find_midpoint(c,d), NE, blue+fontsize(n));
	label("-0.11545", find_midpoint(d,e), NW, blue+fontsize(n));
	label("-0.0816", find_midpoint(e,a), SW, blue+fontsize(n));

	pair ln=(e.x-0.5*size,e.y);
	pair lf=(e.x-0.5*size*2,e.y);
	pair rn=(c.x+0.5*size,c.y);
	pair rf=(c.x+0.5*size*2,c.y);

	//draw(a--ln--d, black+dotted);
	//draw(a--lf--d, black+dotted);
	//draw(b--rn--d, black+dotted);
	//draw(b--rf--d, black+dotted);


	for(int i=0;i<5;++i){
		//filldraw(shift(pent_vert[i])*scale(size_circle)*unitcircle, white, black+linewidth(width_line));
	}
}

void draw_hex(pair ce=(0,0), real size=1.0){

	pair a=(ce.x-0.5*sqrt(3)*size,ce.y-0.5*size);
	pair b=(ce.x,ce.y-size);
	pair c=(ce.x+0.5*sqrt(3)*size,ce.y-0.5*size);
	pair d=(ce.x+0.5*sqrt(3)*size,ce.y+0.5*size);
	pair e=(ce.x,ce.y+size);
	pair f=(ce.x-0.5*sqrt(3)*size,ce.y+0.5*size);
	draw(a--b--c--d--e--f--cycle, black+linewidth(lw));

	pair hex_vert[]={a,b,c,d,e,f};

	//EDIT VERTEX LABELS OF HEXAGON HERE
	label("YBR262C", a, 0.5N+W, black+fontsize(n));
	label("YGR163W", b, S, black+fontsize(n));
	label("YOR221C", c, 0.5N+E, black+fontsize(n));
	label("YOR328W", d, 0.5N+E, black+fontsize(n));
	label("$\textbf{YPL115C}$", e, 4N, black+fontsize(n));
	label("YLR295C", f, 0.5N+W, black+fontsize(n));

	//EDIT VERTEX LABELS OF HEXAGON HERE
	label("1", a, 2S+6W, brown+fontsize(n));
	label("2", b, 4S, brown+fontsize(n));
	label("3", c, 2S+7E, brown+fontsize(n));
	label("4", d, 7E+2S, brown+fontsize(n));
	label("$\textbf{5}$", e, N, brown+fontsize(n));
	label("6", f, 6W+2S, brown+fontsize(n));

	//EDIT EDGE LABELS OF HEXAGON HERE
	label("NA", find_midpoint(a,b), SW, blue+fontsize(n));
	label("0.0149", find_midpoint(b,c), SE, blue+fontsize(n));
	label("-0.1699", find_midpoint(c,d), E, blue+fontsize(n));
	label("NA", find_midpoint(d,e), NE, blue+fontsize(n));
	label("NA", find_midpoint(e,f), NW, blue+fontsize(n));
	label("NA", find_midpoint(a,f), W, blue+fontsize(n));

	pair lnu=(f.x-0.5*size,f.y);
	pair lfu=(f.x-0.5*size*2,f.y);
	pair lnd=(a.x-0.5*size,a.y);
	pair lfd=(a.x-0.5*size*2,a.y);
	pair rnu=(d.x+0.5*size,d.y);
	pair rfu=(d.x+0.5*size*2,d.y);
	pair rnd=(c.x+0.5*size,c.y);
	pair rfd=(c.x+0.5*size*2,c.y);

	//draw(e--lnu--lnd--b, black+dotted);
	//draw(e--lfu--lfd--b, black+dotted);
	//draw(e--rnu--rnd--b, black+dotted);
	//draw(e--rfu--rfd--b, black+dotted);


	for(int i=0;i<6;++i){
		//filldraw(shift(hex_vert[i])*scale(size_circle)*unitcircle, white, black+linewidth(width_line));
	}
}
//---------------Definitions of functions end-------------------------------------------------------------------------------


// DRAWING SHAPES



//triangle 
draw_tri((0.6,1.0), 1.2); //first argument is (x,y) position and second argument is the size
// first argument (0.6,1.2) is the coordinate and second argument 0.75 is the size

//square 
draw_sq((3.3,1.5), 0.9); //first argument is (x,y) position and second argument is the size

//pentagon 
draw_pent((6.72,1.0), 0.7); //first argument is (x,y) position and second argument is the size

//hexagon 
draw_hex((10,1.5), 0.65); //first argument is (x,y) position and second argument is the size


