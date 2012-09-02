#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*** Type definitions ***/

typedef struct{double x,y,z;} t_3d;  /* 3-D POINTS OR VECTORS */

typedef struct{double r;       /* radius */
			 double x,y,z;   /* position */
			}o_sphere;       /* SPHERE PARAMETERS */

typedef struct{int sidehit;          /* side intersected */
			 double xs,ys,zs;      /* size of sides */
			 double x,y,z;         /* center position */
			}o_box;                       /* BOX PARAMETERS */

typedef struct{t_3d nrm;    /* unit normal */
			 double d;           /* plane constant */
			 t_3d e1, e2, e3;    /* edge vectors */
			 double d1, d2, d3;  /* plane constants */
			 t_3d p1,p2,p3;
			}o_triangle;         /* TRIANGLE PARAMETERS */

typedef struct{  t_3d nrm;
					  double d;
			}o_plane;

typedef struct{double cx,cy,cz;   /* center of surface */
			 double A,B,C,D,E,F,G,H,I,J;  /* coeff */
			 double zmin,zmax;
			}o_superq;         /* SUPERQUADRIC PARAMETERS */


typedef struct{int id;               /* object number */
			 int objtyp;           /* object type */
			 int surfnum;          /* surface number */
			 union {
		 o_sphere    *p_sphere;
		 o_box       *p_box;
		 o_triangle  *p_triangle;
		 o_superq    *p_superq;
		 o_plane		 *p_plane;
		     } objpnt;       /* GENERIC OBJECT TYPE PTR */
	       }t_object;

typedef struct{double x,y,z,bright;} t_light;  /* light source */

typedef struct{double ar,ag,ab;  /* ambient r,g,b */
	       double dr,dg,db;  /* diffuse r,g,b */
			 double sr,sg,sb;  /* specular r,g,b */
	       double coef;      /* specular coef */
	       double refl;      /* reflection 0-1 */
			 double transp;    /* transparency 0-1 */
			 double medium;
			 int pattern;
			}t_surface;        /* SURFACE PROPERTIES PARAMETERS */

typedef struct {double r,g,b;} t_color;  /* COLOR VECTOR */


/*** Constants ***/

#define LIGHTS        4    /* Max # of light sources */
#define OBJECTS       50   /* Max # of objects*/
#define SURFACES      50   /* MAx # of surface types */
#define SCREENWIDTH   320  /* Max screen dimensions */
#define SCREENHEIGHT  200
#define ASPECTRATIO   1.6
/*** #define GAMMA         1.8 ***/
#define OTYPSPHERE    0    /* Object type numbers */
#define OTYPBOX       1
#define OTYPTRIANGLE  2
#define OTYPSUPERQ    3
#define OTYPPLANE     4


/*** Function declarations ***/

double brightness(int source,int lnum,t_3d *pos,t_3d *ray);     /* calculate perceived brightness of source */
int calc_color_index(double intensity);  /* scales an an rgb color to a pixel-plane index */
void copy_v(t_3d *a,t_3d *b);
void scaler(double c,t_3d *v);
void sum(t_3d a,t_3d b,t_3d *c);
void diff(t_3d a,t_3d b,t_3d *c);
void crossp(t_3d *a,t_3d *b,t_3d *o);            /* calculate vector cross produce */
double dotp(t_3d *a,t_3d *b);       		/* calculate dot product */
void endpic(void);            /* close raw pixel file calculated by program */
double intersect(int source,t_3d *pos,t_3d *ray,t_color *color);      /* main program's closest-intersection routine */
void lightray(int lnum,t_3d *objpos,t_3d *lray);          /* compute ray vector from int-pt. to light-src. */
void linepic(int lineno,double pixels[SCREENWIDTH][3]);           /* output a line of pixel colors */
void normalize(t_3d *a);		        /* convert a vector to a unit vector */
double magnitude(t_3d a);
void shade(t_3d *pos,t_3d *ray,t_3d *nrm,t_object *obj,t_color *color);   /* use illum-model to calculate color at int-pt */
void startpic(char outfile[60],int y,int x);          /* open raw pixel file */
void viewing(t_3d *scrnx,t_3d *scrny,t_3d *firstray);           /* calculate 1st ray using viewing parameters */

/* object routines: */
void maksph(int surf,double r,double x,double y,double z);       /* make a sphere */
double intsph(t_3d *pos,t_3d *ray,t_object *obj);    /* calculate ray-sphere intersection point */
void nrmsph(t_3d *pos,t_object *obj,t_3d *nrm);       /* calculate sphere surface normal vector */
void makbox(int surf,double x,double y,double z,double xs,double yx,double zs);
double intbox(t_3d *pos,t_3d *ray,t_object *obj);
void nrmbox(t_3d *pos,t_object *obj,t_3d *nrm);
void maktri(int surf,t_3d p1,t_3d p2,t_3d p3);       					/* Same for a triangles */
double inttri(t_3d *pos,t_3d *ray,t_object *obj);
void nrmtri(t_3d *pos,t_object *obj,t_3d *nrm);
void makplane(int surf,t_3d nrm,t_3d pt);
double intplane(t_3d *pos,t_3d *ray,t_object *obj);
void nrmplane(t_3d *pos,t_object *obj,t_3d *nrm);
void patsph(t_3d pt,t_object *obj,double *u,double *v);
void patbox(t_3d pt,t_object *obj,double *u,double *v);
void pattri(t_3d pt,t_object *obj,double *u,double *v);
void patsup(t_3d pt,t_object *obj,double *u,double *v);
void patplane(t_3d pt,t_object *obj,double *u,double *v);
void checkerbrd(double u,double v,int sizex,int sizey,t_color *color);
void maksup(int surf,double cx,double cy,double cz,
	    double A,double B,double C,double D,double E,double F,
	    double G,double H,double I,double J,double zmin,double zmax);
double intsup(t_3d *pos,t_3d *ray,t_object *obj);
void nrmsup(t_3d *pos,t_object *obj,t_3d *nrm);

void setup(void);       /* Sets up the scene */

/* generic object intersection routine--pointer to appropriate function */
double (*objint[5])(t_3d *pos,t_3d *ray,t_object *obj) = {intsph,intbox,inttri,
							intsup,intplane};

/* generic object normal calculation routine--pointer to appropriate function */
void  (*objnrm[5])(t_3d *pos,t_object *obj,t_3d *nrm) = {nrmsph,nrmbox,nrmtri,
							nrmsup,nrmplane};

void (*objpat[5])(t_3d pt,t_object *obj,double *u,double *v) = {patsph,patbox,
							pattri,patsup,patplane};


/* Global variables: */

int nlight;                    /* # of lights presently in use */
int lightlim = LIGHTS;         /* maximum declared */
t_light light[LIGHTS];         /* array of lights */
int nobject=0;                   /* # of objects presently in use */
int objectlim = OBJECTS;       /* maximum declared */
t_object object[OBJECTS];      /* array of objects */
int nsurface;                  /* # of surfaces presently in use */
int surfacelim = SURFACES;     /* maximum declared */
t_surface surface[SURFACES];   /* array of surfaces */
int sizex, sizey;              /* image sizes */
t_3d eyep, lookp, up;          /* view definition */
double hfov, vfov;             /* fields of view */
int level, maxlevel;           /* recursion level */
char outfilename[60];         /* pixel file name */
t_color background;            /* background color */
double max_intensity;          /* should be calculated in pgm, not here */

#define DEGREETORADIAN (M_PI/180.)
#define FAR_AWAY 99.99E+20
#define MAX_DEPTH 5


/*** The main program: ***/

void main()
{
   int line_y, pixel_x;
	t_3d scrnx,scrny,firstray,ray;
   t_color color;
	double dis, line[SCREENWIDTH][3];

	setup();
	viewing(&scrnx,&scrny,&firstray);
	startpic(outfilename, sizey, sizex);

	for (line_y=0; line_y<sizey; line_y++)
	{
	  for (pixel_x=0; pixel_x<sizex; pixel_x++)
     {
       ray.x = firstray.x + pixel_x*scrnx.x - line_y*scrny.x;
		 ray.y = firstray.y + pixel_x*scrnx.y - line_y*scrny.y;
		 ray.z = firstray.z + pixel_x*scrnx.z - line_y*scrny.z;
		 normalize(&ray);

		 /* actual ray trace */
		 dis = intersect(-1,&eyep,&ray,&color);
		 if (dis>0)     /* ray intersected object  */
		 {
	  line[pixel_x][0] = color.r;
	  line[pixel_x][1] = color.g;
	  line[pixel_x][2] = color.b;
		 }
		 else           /* use background color  */
		 {
	 line[pixel_x][0] = background.r;
	 line[pixel_x][1] = background.g;
	 line[pixel_x][2] = background.b;
		 }
	  }
	  linepic(line_y,line);   /*output line of pixels  */
     if (line_y%10 == 0)
	  {
	printf("done line %d\n",line_y);
	fflush(stdout);
	  }
	}
	endpic();   /*  done with picture  */
}

/*************************************************************/
/***   setup.c -- hardcoded viewing/model parameters:      ***/
/***   These should be read in from a scene file; could be ***/
/***   prepared interactively from the project 1 program   ***/
/*************************************************************/

void setup()
{
		t_3d p1,p2,p3;
		t_3d nrm,pt;

/* set viewing/lighting/recursion-depth/output-file parameters */
   level = 0;              /* recursion level--only for recursive RT */
	sizex = 320.;            /* screen dimensions */
	sizey = 200.;
	hfov = 60.;              /* field of view parameters */
	vfov = 50.;
	eyep.x = 100.0;         /* viewpoint position */
	eyep.y = 0.0;
	eyep.z =100.0;
	lookp.x = 0.0;          /* focus point position */
	lookp.y = 0.0;
	lookp.z = 0.0;
	up.x = 0.0;             /* up direction */
	up.y = 0.0;
	up.z = 1.0;
	strcpy(outfilename,"r0.pix");  /* raw pixel file name */
	nlight = 1;             /* number of light sources */
	light[0].x = 100.;       /* light source position */
	light[0].y = 150.;
	light[0].z = 50.;
	light[0].bright = 5.0;  /* light source intensity */
/* Object definitions follow (parameters could be read from scene file) */

	/* sphere -- surface type 0, radius 40, */
	/* center at (0,0,0) */
	maksph(0,40.0,0.0,0.0,0.0);   /* make the sphere from input parameters */
	/* sphere -- surface type 1 radius 40, */
	/* center at (-45,30,30) */
	maksph(1,40.0,-45.0,30.0,30.0);
	maksph(2,25.,-50.,-70.,60.);

	p1.x=0.;    p1.y=-80.;  p1.z=10.;
	p2.x=-5.;     p2.y=-40.;  p2.z=15.;
	p3.x=-10.;    p3.y=-60.;  p3.z=70.;

	makbox(3,-10.,70.,50.,40.,40.,40.);
	maktri(4,p1,p2,p3);

	maksup(5,20.0,50.0,0.0,1./100.,0.,0.,0.,1./100.,0.,0.,-1./(800.),0.,0.,
			 0.0,70.0);

/* surface properties follow--these should be read in from the scene file */

	nsurface = 6;         /* number of surfaces */

	surface[0].ar = 100.;  /* Surface 0 -- ambient reflection coefficients */
	surface[0].ag = 0.;  /* has maximum red, no green, no blue */
	surface[0].ab = 0.;
	surface[0].dr = 100.;  /* surface 0 -- diffuse reflection coefficients */
	surface[0].dg = 0.;    /* has maximum red, no green, no blue */
	surface[0].db = 0.;
	surface[0].sr = 190.;  /* surface 0 -- specular reflection coefficients */
	surface[0].sg = 190.;  /* has equal mix of red, green, blue */
	surface[0].sb = 190.;
	surface[0].coef = 30.; /* surface 0 specular exponent--highly reflective */
	surface[0].refl=0.;
	surface[0].transp=0.8;
	surface[0].medium=1.0;
	surface[0].pattern=0;

	surface[1].ar = 0.;    // surface 1 -- green duller
	surface[1].ag = 100.;
	surface[1].ab = 0.;
	surface[1].dr = 0.;
	surface[1].dg = 100.;
	surface[1].db = 0.;
	surface[1].sr = 90.;   // specular reflection quite small
	surface[1].sg = 90.;
	surface[1].sb = 90.;
	surface[1].coef = 25.;
	surface[1].refl=0.;
	surface[1].transp=0.;
	surface[1].pattern=0;
	surface[1].medium=1.;

	surface[2].ar = 30.;    // surface 3 - white
	surface[2].ag = 30.;
	surface[2].ab = 30.;
	surface[2].dr = 30.;
	surface[2].dg = 30.;
	surface[2].db = 30.;
	surface[2].sr = 100.;   // specular reflection quite large
	surface[2].sg = 100.;
	surface[2].sb = 100.;
	surface[2].coef = 40.;
	surface[2].refl=0.9;
	surface[2].transp=0.;
	surface[2].pattern=0;
	surface[2].medium=1.0;

	surface[3].ar =0.;
	surface[3].ag=0.;
	surface[3].ab=100.;
	surface[3].dr=0.;
	surface[3].dg=0.;
	surface[3].db=100.;
	surface[3].sr=200.;
	surface[3].sg=200.;
	surface[3].sb=200.;
	surface[3].coef = 30.;
	surface[3].refl = 0.;
	surface[3].transp=0.;
	surface[3].medium=1.0;
	surface[3].pattern=1.;

	surface[4].ar=300.;
	surface[4].ag=20.;
	surface[4].ab=0.;
	surface[4].dr=300.;
	surface[4].dg=20.;
	surface[4].db=0.;
	surface[4].sr=100.;
	surface[4].sg=100.;
	surface[4].sb=100.;
	surface[4].coef=30.;
	surface[4].refl=0.;
	surface[4].transp=0.2;
	surface[4].pattern=1;
	surface[4].medium=1.0;

	surface[5].ar=150.;
	surface[5].ag=0.;
	surface[5].ab=0.;
	surface[5].dr=150.;
	surface[5].dg=0.;
	surface[5].db=0.;
	surface[5].sr=200.;
	surface[5].sg=200.;
	surface[5].sb=200.;
	surface[5].coef=25.;
	surface[5].refl=0.;
	surface[5].transp=0.5;
	surface[5].medium=1.5;
	surface[5].pattern=0;

	background.r = 250.;   /* background color */
	background.g = 250.;   /* gray   */
	background.b = 250.;

	max_intensity = 390.;  /* This should be calculated by the program */
}

/*** viewing function--calculate screen coordinate axis vectors    ***/
/*** scrnx & scrny), scaled to width and height of a pixel, and    ***/
/*** the ray vector from eyepoint to the top-left pixel (firstray) ***/

double radian(double degree)
{
	return degree*M_PI/180.;
}

void viewing(t_3d *scrnx,t_3d *scrny,t_3d *firstray)
{
	t_3d G;
	t_3d xp,yp,fr;
	double mag_x,mag_y;
	diff(lookp,eyep,&G);
	crossp(&G,&up,scrnx);
	copy_v(scrnx,&xp);
	crossp(&xp,&G,scrny);
	normalize(scrnx);	normalize(scrny);
	mag_y=2.*magnitude(G)*tan(radian(vfov))/200.;
	mag_x=2.*magnitude(G)*tan(radian(hfov))/320.;
	scaler(mag_y,scrny);
	scaler(mag_x,scrnx);
	copy_v(scrnx,&xp);   copy_v(scrny,&yp);
	scaler(160.,&xp);
	scaler(100.,&yp);
	diff(G,xp,firstray);
	copy_v(firstray,&fr);
	sum(fr,yp,firstray);
}

/**** The following functions relate to file output ****/

int width;
FILE *filept;

void startpic(char outfile[60],int y,int x)
{
   width = x;
	filept = fopen(outfile,"wb");
   fwrite(&x,sizeof(int),1,filept);
	fwrite(&y,sizeof(int),1,filept);
}

void linepic(int lineno,double pixels[SCREENWIDTH][3])
{
   unsigned char buffer[3][SCREENWIDTH];
   int i,r,g,b;
	double dr,dg,db;

	for (i=0; i<width; i++)
   {
		r = calc_color_index(pixels[i][0]);
		g = calc_color_index(pixels[i][1]);
      b = calc_color_index(pixels[i][2]);
      buffer[0][i] = r;
      buffer[1][i] = g;
      buffer[2][i] = b;
	}
   fwrite(&lineno,sizeof(int),1,filept);
   fwrite(buffer[0],sizeof(char),width,filept); /* red pixels for a line */
	fwrite(buffer[1],sizeof(char),width,filept); /* green pixels for a line */
   fwrite(buffer[2],sizeof(char),width,filept); /* blue pixels for a line */
}

void endpic(void)
{
   fclose(filept);
}

/* calc_color_intensity()--scale an rgb intensity to 5 bits (32 shades) */

int calc_color_index(double intensity)
{
	int ival;
	double dval;

	/* scale to 0-1 range */
	dval = intensity/max_intensity;   /* max_intensity should be computed */
	if (dval > 1.0) dval = 1.0;
	if (dval < 0.0) dval = 0.0;
	/* convert to integer 0-31 --- 5 bits per channel hardware color mode */
	dval *= 31.0;
	ival = (int) (dval + 0.5);
   return(ival);
}


/* generic object intersection routine--Finds closest intersection point  */
/* between a ray & any object in scene. Input parameters: source object   */
/* ID #, pointers to: ray's starting point and the ray vector. Calculates */
/* position of closest intersection point, normal vector at that point,   */
/* and calls shade() to compute the color of the point. Returns distance  */
/* to intersection point                                                  */

double intersect(int source,t_3d *pos,t_3d *ray,t_color* color)
{
	int objhit,objtry;
   double s,ss;
	t_3d hit,normal;
	t_surface *surf;
	double u,v;

	objhit = -1;
	ss = FAR_AWAY;
	/* check for intersection of ray with all objects */
	for (objtry=0; objtry<nobject; objtry++)
	{
		if (objtry != source)  /* don't try source */
		{
	 s = (*objint[object[objtry].objtyp])
			 (pos,ray,&object[objtry]);
	 /* keep track of closest intersection */
	 if ((s > 0.0) && (s < ss))
	 {
		 objhit = objtry;
		 ss = s;
	 }
		}
	}
	if (objhit < 0) return (0);  /* ray hit no objects */

   /* find point of intersection */
	hit.x = pos->x + ss * ray->x;
   hit.y = pos->y + ss * ray->y;
   hit.z = pos->z + ss * ray->z;
	/* find normal */
	(*objnrm[object[objhit].objtyp])
	  (&hit,&object[objhit],&normal);
	/* find color at point of intersection */
	shade(&hit,ray,&normal,&object[objhit],color);
	surf = &surface[object[objhit].surfnum];
	if (surf->pattern == 1){
		(*objpat[object[objhit].objtyp])
			(hit,&object[objhit],&u,&v);
		checkerbrd(u,v,10,10,color);
	}
	return(ss);
}
/* non-recursive shading routine -- uses the Phong illumination model   */
/* to calculate the color of an intersection point between a ray and an */
/* object. Input parameters: pointers to: position of intersection      */
/* point, the "incident" ray vector, the normal vector at the           */
/* intersection point, and the intersecting object structure. Result    */
/* is a pointer to an rgb color. This will have to be modified for the  */
/* recursive raytracer, taking into account reflection & transparency.  */

void shade(t_3d *pos,t_3d *ray,t_3d *nrm,t_object *obj,t_color *color)
{
	int lnum;
	double k,bright,spec,diffuse;
	double nr,alpha,beta,cos_t;
	t_surface *surf;
	t_3d refl,trans,ltray,I;
	t_color newcolor;
	double dist;
	int id=obj->id;

	/* ambient light contribution */
	surf = &surface[obj->surfnum];
	color->r = surf->ar;
	color->g = surf->ag;
	color->b = surf->ab;

	/* calculate reflected ray */
	k = -2.0 * dotp(ray,nrm);
	refl.x = k*nrm->x + ray->x;
	refl.y = k*nrm->y + ray->y;
	refl.z = k*nrm->z + ray->z;
	if ((surf->refl > 0) && (level <= MAX_DEPTH))
	{
		level++;
		dist=intersect(id,pos,&refl,&newcolor);
		if (dist > 0){
			color->r+=(newcolor.r*surf->refl);
			color->g+=(newcolor.g*surf->refl);
			color->b+=(newcolor.b*surf->refl);
		}
		else {
			color->r+=(background.r*surf->refl);
			color->g+=(background.g*surf->refl);
			color->b+=(background.b*surf->refl);
		}
		level--;
	}

	if ((surf->transp > 0) && (level <=MAX_DEPTH ))
	{
		level++;
		nr=1./surf->medium;
		if (nr>1.0) nr=1.0;
		alpha=nr;
		cos_t=-dotp(ray,nrm);
		beta=nr*cos_t-sqrt(1.-nr*nr*(1.-cos_t*cos_t));
		trans.x=alpha*ray->x+beta*nrm->x;
		trans.y=alpha*ray->y+beta*nrm->y;
		trans.z=alpha*ray->z+beta*nrm->z;
		dist=intersect(id,pos,&trans,&newcolor);
		if (dist > 0){
		   color->r=color->r*(1.-surf->transp)+(newcolor.r*surf->transp);
		   color->g=color->g*(1.-surf->transp)+(newcolor.g*surf->transp);
		   color->b=color->b*(1.-surf->transp)+(newcolor.b*surf->transp);
		}
		else{
		   color->r=color->r*(1.-surf->transp)+(background.r*surf->transp);
		   color->g=color->g*(1.-surf->transp)+(background.g*surf->transp);
		   color->b=color->b*(1.-surf->transp)+(background.b*surf->transp);
		}
		level--;
	}

	for (lnum=0; lnum<nlight; lnum++)
	{
		/* get ray to light */
      lightray(lnum,pos,&ltray);
		diffuse = dotp(nrm,&ltray);
      if (diffuse > 0)
      {
	 /* object faces light, add diffuse */
	 bright = brightness(obj->id,lnum,pos,&ltray);
	 diffuse *= bright;
	 color->r += surf->dr * diffuse;
	 color->g += surf->dg * diffuse;
	 color->b += surf->db * diffuse;

	 spec = dotp(&refl,&ltray);
	 if (spec > 0)
	 {
	    /* highlight is here, add specular */
	    spec = bright * pow(spec,surf->coef);
	    color->r += surf->sr * spec;
	    color->g += surf->sg * spec;
	    color->b += surf->sb * spec;
	 }
		}
	}
}

/*** lighting functions lightray() & brightness() ***/

/* lightray() calculates vector from object point objpos to light source */
/* number lnum                                                           */

void lightray(int lnum,t_3d *objpos,t_3d *lray)
{
	lray->x = light[lnum].x - objpos->x;
	lray->y = light[lnum].y - objpos->y;
	lray->z = light[lnum].z - objpos->z;
	normalize(lray);
}

/* brightness() returns intensity of light source # lnum as seen from  */
/* position pos along ray direction -- Trivial for simple tracer, but  */
/* is a shadow feeler for the recursive ray tracer. Will need to make  */
/* calls to object intersection routines. If ray is blocked, returned  */
/* brightness is zero.                                                 */

double brightness(int source,int lnum,t_3d *pos,t_3d *ray)
{
	int i;
	float t;
	float ss,s;
	t_3d dist_betw;
	t = (light[lnum].x - pos->x)/ray->x;
	if (t==0)
		t = (light[lnum].y - pos->y)/ray->y;
	if (t==0)
		t = (light[lnum].z - pos->z)/ray->z;

	ss = FAR_AWAY;
	/* check for intersection of ray with all objects */
	for (i=0; i<nobject; i++)
	{
		if (i != source)  /* don't try source */
		{
		s = (*objint[object[i].objtyp])
			  (pos,ray,&object[i]);
	 /* keep track of closest intersection */
	 if ((s > 0.0) && (s < ss))
		 ss = s;
		}
	}
	if (ss < t) return 0.0;
	return(light[lnum].bright);  
}

void checkerbrd(double u,double v,int sizex,int sizey,t_color *color)
{
	int T;
	double brighter=100.;
	T=(int)((float)sizex*u+0.5) + (int)((float)sizey*v+0.5);
	if (T%2==0){
		color->r +=brighter;
		color->g += brighter;
		color->b +=brighter;
	}
	else{
		color->r -=brighter;
		color->g -= brighter;
		color->b -=brighter;
	}
}

/*** vector math utility routines ***/
void copy_v(t_3d *a,t_3d* b)
{
	b->x=a->x;
	b->y=a->y;
	b->z=a->z;
}

void normalize(t_3d *a)
{
	double d;
	d = sqrt ( a->x * a->x + a->y * a->y + a->z * a->z );
	a->x /= d;
	a->y /= d;
	a->z /= d;
}

double magnitude(t_3d a)
{
	return sqrt ( a.x * a.x + a.y * a.y + a.z * a.z );
}

void scaler(double c,t_3d* v)
{
	v->x*=c;
	v->y*=c;
	v->z*=c;
}

void diff(t_3d a,t_3d b,t_3d* c)
{
	c->x=a.x-b.x;
	c->y=a.y-b.y;
	c->z=a.z-b.z;
}

void sum(t_3d a,t_3d b,t_3d* c)
{
	c->x=a.x+b.x;
	c->y=a.y+b.y;
	c->z=a.z+b.z;
}

double dotp(t_3d *a,t_3d *b)
{
	double d;
	d = ( a->x * b->x + a->y * b->y + a->z * b->z );
	return(d);
}

void crossp(t_3d *a,t_3d *b,t_3d *o)
{
	o->x = (a->y * b->z) - (a->z * b->y);
	o->y = (a->z * b->x) - (a->x * b->z);
	o->z = (a->x * b->y) - (a->y * b->x);
}


/*** sphere functions ***/

/* create sphere object */

void maksph(int surf,double r,double x,double y,double z)
{
   o_sphere *sphere;

	sphere = (o_sphere *) malloc(sizeof(o_sphere));
   object[nobject].id = nobject;
	object[nobject].objtyp = OTYPSPHERE;
   object[nobject].surfnum = surf;
   object[nobject].objpnt.p_sphere = sphere;
   sphere->r = r;
   sphere->x = x;
   sphere->y = y;
   sphere->z = z;
   nobject++;
}

double squared(double a)
{
  return a*a;
}


/* find roots of a quadratic polynomial */
int solve_quad(double a,double b,double c,double *r1,double *r2)
{
   double discrim=squared(b)-4.*a*c;

   if (a==0) return -1;
   if (discrim < 0.) return -1;

	*r1=(-b+sqrt(discrim))/(2.*a);
	*r2=(-b-sqrt(discrim))/(2.*a);
	return 1;
}

double min_pos(double n1,double n2)
{
   if ((n1 <=0 ) && (n2 <=0)) return 0;
   if (n1<=0) return n2;
   if (n2 <=0) return n1;
   if (n1 < n2) return n1;
   return n2;
}

/* calculate ray-sphere intersection point */
double intsph(t_3d *pos,t_3d *ray,t_object *obj)
{
   o_sphere *sphere;
   t_3d center;
   double radius;
   double C1,C2,C3;
   double A,B,C;
	double root1=0,root2=0;

	sphere=obj->objpnt.p_sphere;
   center.x=sphere->x;
   center.y=sphere->y;
   center.z=sphere->z;
   radius=sphere->r;

   C1=squared(radius)-squared(center.x)-squared(center.y)-squared(center.z);
	C2=C1-squared(pos->x)-squared(pos->y)-squared(pos->z);

	C3=C2+2.*pos->x*center.x+2.*pos->y*center.y+
                           2.*pos->z*center.z;

	A=squared(ray->x)+squared(ray->y)+squared(ray->z);
	B=2.*(pos->x*ray->x-ray->x*center.x+pos->y*ray->y-ray->y*center.y+
			pos->z*ray->z-ray->z*center.z);
   C=-C3;
   if (!solve_quad(A,B,C,&root1,&root2))
	return FAR_AWAY;
   return min_pos(root1,root2);
}

/* calculate unit normal to sphere at intersection point */

void nrmsph(t_3d *pos,t_object *obj,t_3d *nrm)
 /*  t_3d *pos;                  point of intersection
	  t_object *obj;              sphere description
	  t_3d *nrm;                  return surface normal */

{
  o_sphere *sphere;
  sphere=obj->objpnt.p_sphere;
  nrm->x=pos->x-sphere->x;
  nrm->y=pos->y-sphere->y;
  nrm->z=pos->z-sphere->z;
  normalize(nrm);
}

/*** The following object routine stubs need to be implemented ***/

void makbox(int surf,double x,double y,double z,double xs,double ys,double zs)
{
    o_box *box;

    box=(o_box *)malloc(sizeof(o_box));
    object[nobject].id=nobject;
    object[nobject].objtyp=OTYPBOX;
    object[nobject].surfnum=surf;
    object[nobject].objpnt.p_box=box;
    box->x=x;	box->y=y;	box->z=z;
    box->xs=xs;	box->ys=ys;	box->zs=zs;
    box->sidehit=-1;
	 nobject++;
}

void makplane(int surf,t_3d nrm,t_3d pt)
{
	o_plane *plane;
	plane=(o_plane *)malloc(sizeof(o_plane));
	object[nobject].id=nobject;
	object[nobject].objtyp=OTYPPLANE;
	object[nobject].surfnum=surf;
	object[nobject].objpnt.p_plane=plane;

	normalize(&nrm);
	plane->nrm.x=nrm.x; plane->nrm.y=nrm.y;  plane->nrm.z=nrm.z;
	plane->d=-dotp(&nrm,&pt);
	nobject++;
}
void maktri(int surf,t_3d p1,t_3d p2,t_3d p3)
{
	t_3d v1,v2,v3;
	o_triangle *triangle;

	triangle=(o_triangle *)malloc(sizeof(o_triangle));
	object[nobject].id=nobject;
	object[nobject].objtyp=OTYPTRIANGLE;
	object[nobject].surfnum=surf;
	object[nobject].objpnt.p_triangle=triangle;

	copy_v(&p1,&(triangle->p1));
	copy_v(&p2,&(triangle->p2));
	copy_v(&p3,&(triangle->p3));
	diff(p2,p1,&v1);
	diff(p3,p2,&v2);
	diff(p1,p3,&v3);

	crossp(&v1,&v2,&(triangle->nrm));
	normalize(&(triangle->nrm));
	crossp(&(triangle->nrm),&v1,&(triangle->e1));
	crossp(&(triangle->nrm),&v2,&(triangle->e2));
	crossp(&(triangle->nrm),&v3,&(triangle->e3));
	normalize(&(triangle->e1));
	normalize(&(triangle->e2));
	normalize(&(triangle->e3));
	triangle->d=-dotp(&(triangle->nrm),&p1);
	triangle->d1=dotp(&p1,&(triangle->e1));
	triangle->d2=dotp(&p2,&(triangle->e2));
	triangle->d3=dotp(&p3,&(triangle->e3));
	nobject++;
}

void maksup(int surf,double cx,double cy,double cz,double A,double B,
		double C,double D,double E,double F,double G,double H,
		double I,double J,double zmin,double zmax)
{
	o_superq *superq;
	superq=(o_superq *)malloc(sizeof(o_superq));
	object[nobject].id=nobject;
	object[nobject].objtyp=OTYPSUPERQ;
	object[nobject].surfnum=surf;
	object[nobject].objpnt.p_superq=superq;

	superq->cx=cx; superq->cy=cy; superq->cz=cz;
	superq->A=A; superq->B=B; superq->C=C;
	superq->D=D; superq->E=E; superq->F=F;
	superq->G=G; superq->H=H; superq->I=I;
	superq->J=J; superq->zmin=zmin; superq->zmax=zmax;
	nobject++;
}

int between(double b,double a,double c)
{
	if ((b<=c)&&(b>=a)) return 1;
	return 0;
}

double intbox(t_3d *pos,t_3d *ray,t_object *obj)
{
	o_box *box;
	double t[6];
	double xc,yc,zc;
	double xs,ys,zs;
	double min_val;
	int min;
	int i;
	t_3d point[6];

	box=obj->objpnt.p_box;
	obj->objpnt.p_box->sidehit=-1;
	xc=box->x;  yc=box->y; zc=box->z;
	xs=box->xs; ys=box->ys; zs=box->zs;

	if (ray->z == 0.0){
	   t[0]=FAR_AWAY; t[5]=FAR_AWAY;
	}
	else{
	   t[0]=(zc-zs/2.0-pos->z)/(ray->z);
		t[5]=(zc+zs/2.0-pos->z)/(ray->z);
	}

	if (ray->x == 0.0){
	  t[1]=FAR_AWAY; t[3]=FAR_AWAY;
	}
	else{
	   t[1]=(xc+xs/2.0-pos->x)/(ray->x);
	   t[3]=(xc-xs/2.0-pos->x)/(ray->x);
	}

	if (ray->y == 0.0){
	  t[2]=FAR_AWAY; t[4]=FAR_AWAY;
	}
	else{
	   t[2]=(yc+ys/2.0-pos->y)/(ray->y);
	   t[4]=(yc-ys/2.0-pos->y)/(ray->y);
	}

	for (i=0;i<6;i++){
	   point[i].x=pos->x+t[i]*ray->x;
	   point[i].y=pos->y+t[i]*ray->y;
	   point[i].z=pos->z+t[i]*ray->z;
	}

	if (!between(point[0].x,xc-xs/2.0,xc+xs/2.0) ||
		  !between(point[0].y,yc-ys/2.0,yc+ys/2.0)){
		t[0]=FAR_AWAY;
	     }
	if (!between(point[1].y,yc-ys/2.0,yc+ys/2.0) ||
	    !between(point[1].z,zc-zs/2.0,zc+zs/2.0)){
		t[1]=FAR_AWAY;
	    }
	if (!between(point[2].x,xc-xs/2.0,xc+xs/2.0) ||
	    !between(point[2].z,zc-zs/2.0,zc+zs/2.0)){
		t[2]=FAR_AWAY;
	    }
	if (!between(point[3].y,yc-ys/2.0,yc+ys/2.0) ||
	    !between(point[3].z,zc-zs/2.0,zc+zs/2.0)){
		t[3]=FAR_AWAY;
		 }
	if (!between(point[4].x,xc-xs/2.0,xc+xs/2.0) ||
	    !between(point[4].z,zc-zs/2.0,zc+zs/2.0)){
		t[4]=FAR_AWAY;
	    }
	if (!between(point[5].x,xc-xs/2.0,xc+xs/2.0) ||
	    !between(point[5].y,yc-ys/2.0,yc+ys/2.0)){
		t[5]=FAR_AWAY;
	    }
	for (i=0;i<6;i++){
	   if (t[i]<=0.0) t[i]=FAR_AWAY;
	}

		 min_val=FAR_AWAY;  min=-1;
       for (i=0;i<6;i++){
	  if (t[i]<=min_val) {min_val=t[i]; min=i;}
       }
       if (min_val<FAR_AWAY){
	  obj->objpnt.p_box->sidehit=min;
	  return min_val;
       }
       return FAR_AWAY;
}

int inside(t_3d point,o_triangle *triangle)
{
	if (dotp(&point,&(triangle->e1)) < triangle->d1) return 0;
	if (dotp(&point,&(triangle->e2)) < triangle->d2) return 0;
	if (dotp(&point,&(triangle->e3)) < triangle->d3) return 0;
	return 1;
}

double intplane(t_3d *pos,t_3d *ray,t_object *obj)
{
	o_plane *plane;
	double bottom,top;
	double t;
	plane=obj->objpnt.p_plane;
	bottom=dotp(&(plane->nrm),ray);
	if (bottom==0.) return FAR_AWAY;
	top=-(dotp(&(plane->nrm),pos)+plane->d);
	t=top/bottom;
	if (t<=0.) return 0.0;
	return t;
}

double inttri(t_3d *pos,t_3d *ray,t_object *obj)
{
	o_triangle *triangle;
	double bottom,top;
	double t;
	t_3d point;
	triangle=obj->objpnt.p_triangle;
	bottom=dotp(&(triangle->nrm),ray);
	if (bottom==0.) return FAR_AWAY;
	top=-1.*(dotp(&(triangle->nrm),pos)+triangle->d);
	t=top/bottom;
	if (t<=0.) return 0.0;
	point.x=pos->x+t*(ray->x);
	point.y=pos->y+t*(ray->y);
	point.z=pos->z+t*(ray->z);
	if (!(inside(point,triangle))) return FAR_AWAY;
	return t;
}

double intsup(t_3d *pos,t_3d *ray,t_object *obj)
{
	double K1,K2,K3,K4,K5,K6,K7;
	o_superq *superq;
	double A,B,C,D,E,F,G,H,I,J;
	double cx,cy,cz;
	double x0,y0,z0;
	double xd,yd,zd;
	double root1=0.0,root2=0.0;
	double h1,h2;

	superq=obj->objpnt.p_superq;
	x0=pos->x;   y0=pos->y;  z0=pos->z;
	xd=ray->x;   yd=ray->y;  zd=ray->z;

	cx=superq->cx;  cy=superq->cy; cz=superq->cz;
	A=superq->A; B=superq->B; C=superq->C;
	D=superq->D; E=superq->E; F=superq->F;
	G=superq->G; H=superq->H; I=superq->I;
	J=superq->J;

	K1=-A*cx*cx-2.*B*cx*cy-2.*C*cx*cz+2.*D*cx-E*cy*cy-2.*F*cy*cz-
	   H*cz*cz+2.*I*cz-J;
	K2=2.*A*cx+2.*B*cy+2.*C*cz-2.*D;
	K3=2.*B*cx+2.*E*cy+2.*F*cz-2.*G;
	K4=2.*C*cx+2.*F*cy+2.*H*cz-2.*I;

	K5=K1-A*x0*x0+K2*x0-E*y0*y0+K3*y0-H*z0*z0+K4*z0-2.*B*x0*y0-
	   2.*C*x0*z0-2.*F*y0*z0;
	K6=A*xd*xd+E*yd*yd+H*zd*zd+2.*B*xd*yd+2.*C*xd*zd+2.*F*yd*zd;
	K7=2.*A*x0*xd-K2*xd+2*E*y0*yd-K3*yd+2.*H*z0*zd-K4*zd+2.*B*x0*yd+
	   2.*B*y0*xd+2.*C*x0*zd+2.*C*z0*zd+2.*F*y0*zd+2.*F*z0*yd;
	if (!solve_quad(K6,K7,-K5,&root1,&root2))
		return FAR_AWAY;
	h1=pos->z+root1*ray->z;	h2=pos->z+root2*ray->z;
	if ((h1>superq->zmax) || (h1<superq->zmin)) root1=0.0;
	if ((h2>superq->zmax) || (h2<superq->zmin)) root2=0.0;
	return min_pos(root1,root2);
}

void nrmbox(t_3d *pos,t_object *obj,t_3d *nrm)
{
	o_box *box;
	box=obj->objpnt.p_box;
	switch (box->sidehit){
	 case 0: {
		   nrm->x=0.0; nrm->y=0.0; nrm->z=-1.0;
		   break;
		}
	 case 1: {
		 nrm->x=1.0; nrm->y=0.0; nrm->z=0.0;
		 break;
		}
	 case 2:{
		 nrm->x=0; nrm->y=1.0; nrm->z=0.0;
		 break;
		}
	 case 3:{
		 nrm->x=-1.0; nrm->y=0.0; nrm->z=0.0;
		 break;
		}
	 case 4:{
		 nrm->x=0.0; nrm->y=-1.0; nrm->z=0.0;
		 break;
		}
	 case 5:{
		 nrm->x=0.0; nrm->y=0.0; nrm->z=1.0;
		 break;
		}
	 case -1:{
		printf("bug");
		break;
	 }
	}
}

void nrmtri(t_3d *pos,t_object *obj,t_3d *norm)
{
	o_triangle *triangle;
	triangle=obj->objpnt.p_triangle;
	norm->x=triangle->nrm.x;
	norm->y=triangle->nrm.y;
	norm->z=triangle->nrm.z;
}

void nrmplane(t_3d *pos,t_object *obj,t_3d *nrm)
{
	o_plane *plane;
	plane=obj->objpnt.p_plane;
	nrm->x=plane->nrm.x;
	nrm->y=plane->nrm.y;
	nrm->z=plane->nrm.z;
}

void nrmsup(t_3d *pos,t_object *obj,t_3d *nrm)
{
	o_superq *superq;
	double x,y,z;
	superq=obj->objpnt.p_superq;
	x=pos->x-superq->cx;
	y=pos->y-superq->cy;
	z=pos->z-superq->cz;
	nrm->x=2.*superq->A*x+2.*superq->B*y+2.*superq->C*z+2.*superq->D;
	nrm->y=2.*superq->B*x+2.*superq->E*y+2.*superq->F*z+2.*superq->G;
	nrm->z=2.*superq->C*x+2.*superq->F*y+2.*superq->H*z+2.*superq->I;
	normalize(nrm);
}
void patsph(t_3d pt,t_object *obj,double *u,double*v)
{
	o_sphere *sphere;
	double theta,phi;
	double R,X,Y,Z;
	double xc,yc,zc;
	sphere=obj->objpnt.p_sphere;
	R=sphere->r;
	xc=sphere->x; yc=sphere->y; zc=sphere->z;
	X=(pt.x-xc)/R;  Y=(pt.y-yc)/R; Z=(pt.z-zc)/R;
	if (Z>=1.0) Z=1.0;
	if (Z<=-1.0) Z=-1.0;
	phi=acos(Z);
	*v=(M_PI-phi)/M_PI;
	if (X/sin(phi)>=1.0)
		theta=acos(1.0);
	else if (X/sin(phi)<=-1.0) theta=acos(-1.0);
	else theta=acos(X/sin(phi));
	if (Y > 0.) *u=theta/(2.*M_PI);
	else *u=1.-(theta/(2.*M_PI));
}

void patbox(t_3d pt,t_object *obj,double *u,double*v)
{
	o_box *box;
	double xc,yc,zc;
	double xs,ys,zs;
	double x,y,z;
	double w,h;
	double xp,yp;
	box=obj->objpnt.p_box;

	x=pt.x-box->x;	y=pt.y-box->y; z=pt.z-box->z;
	xs=box->xs;	ys=box->ys; zs=box->zs;
	w=2.*xs+2.*ys;
	h=2.*xs+zs;
	x-=xs/2.0;
	y+=ys/2.0;
	z+=zs/2.0;

	switch (box->sidehit){
		case 0:{
			xp=y;
			yp=xs+x;
			break;
		}
		case 1:{
			xp=y;
			yp=xs+z;
			break;
		}
		case 2:{
			xp=ys-x;
			yp=xs+z;
			break;
		}
		case 3:{
			xp=ys+xs+ys-y;
			yp=xs+z;
			break;
		}
		case 4:{
			xp=ys+xs+ys+xs+x;
			yp=xs+z;
			break;
		}
		case 5:{
			xp=y;
			yp=xs+zs+x;
			break;
		}
	}
	*u=xp/w; *v=yp/h;
}

double dist(t_3d p1,t_3d p2)
{
	return sqrt(squared(p1.x-p2.x)+squared(p1.y-p2.y)+
					squared(p1.z-p2.z));
}

void pattri(t_3d pt,t_object *obj,double *u,double*v)
{
	o_triangle *triangle;
	double x,y,z;
	double w1,w2;
	double k;
	t_3d S,T,V;

	triangle=obj->objpnt.p_triangle;

	diff(triangle->p2,triangle->p1,&S);
	normalize(&S);
	crossp(&(triangle->nrm),&(S),&T);
	V.x=pt.x-triangle->p1.x;
	V.y=pt.y-triangle->p1.y;
	V.z=pt.z-triangle->p1.z;

	w1=dist(triangle->p1,triangle->p2);
	w2=dist(triangle->p1,triangle->p3);
	if (w1>w2) k=w1;
	else k=w2;

	*u=dotp(&V,&S)/k;
	*v=dotp(&V,&T)/k;


}

void patsup(t_3d pt,t_object *obj,double *u,double*v)
{
}

void patplane(t_3d pt,t_object *obj,double *u,double*v)
{
}