settings.outformat="pdf";
texpreamble("\usepackage[charter]{mathdesign}");
//defaultpen(fontsize(1pt));

unitsize(4cm);
dotfactor=10;

real radius=0.075;
real width_circle_circumference =2;
real width_line=1.5;
real height_of_turn=0.5;
real arrow_size=5;


//Coordinates of circles
pair position_circle[] = {(2,0), (4,0)};

string title[]={"X", "Y"};

draw((position_circle[0].x,position_circle[0].y+radius)::(3,height_of_turn)::(position_circle[1].x,position_circle[1].y+radius), black+linewidth(width_line), bar=Bar(arrow_size));
draw((position_circle[1].x,position_circle[1].y-radius)::(3,-height_of_turn)::(position_circle[0].x,position_circle[0].y-radius), black+linewidth(width_line), bar=Bar(arrow_size));

for(int i=0; i<2; ++i){
	//filldraw(shift(position_circle[i])*scale(radius)*unitcircle, white, black+linewidth(width_line));
	label(title[i], position_circle[i]);
}


//clip(shift(0.6,0.6)*scale(2.3)*unitsquare);

//label("$t_{pd}$",(2,2.85),W, cpen);
//label("$t_{pp'}$",(2.1,2.85),E, dpen);
//label("$t_{pp}$",(2.4,2.85),E, open);
