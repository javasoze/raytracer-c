package raytracerj;

import java.util.List;

public class Renderer {

  public static final int DEFAULT_DEPTH = 5;
  
  private int screenWidth;
  
  private int screenHeight;
  private int depth = DEFAULT_DEPTH;
  private Scene scene;

  private Point3d scrnx = new Point3d();
  private Point3d scrny = new Point3d();
  private Point3d firstray = new Point3d();
  
  private int level;
  
  
  public Renderer(int screenWidth,int screenHeight,Scene scene){
    this.screenWidth = screenWidth;
    this.screenHeight = screenHeight;
    this.scene = scene;
  }
  
  public void initScene(){
    double magX,magY;
    Point3d fr = new Point3d();
    Point3d xp = new Point3d();
    Point3d yp = new Point3d();
    Point3d G = new Point3d();
    Point3d.diff(scene.look, scene.eye, G);
    Point3d.crossp(G, scene.up, scrnx);
    Point3d.copy(scrnx, xp);
    Point3d.crossp(xp, G, scrny);
    scrnx.normalize();
    scrny.normalize();
    magY = 2.0 * G.magnitude() * Math.tan(MathUtil.radian(scene.vfov)) / (double)this.screenHeight;
    magX = 2.0 * G.magnitude() * Math.tan(MathUtil.radian(scene.hfov)) / (double)this.screenWidth;
    scrny.scale(magY);
    scrnx.scale(magX);
    Point3d.copy(scrnx, xp);
    Point3d.copy(scrny, yp);
    xp.scale((double)this.screenWidth/2.0);
    yp.scale((double)this.screenHeight/2.0);
    Point3d.diff(G, xp, firstray);
    Point3d.copy(firstray, fr);
    Point3d.sum(fr, yp, firstray);
    level = 0;
  }
  
  public void render(){
    Point3d scrny = new Point3d();
    Point3d firstray = new Point3d();
    Point3d ray = new Point3d();
    Color color = new Color();
    
    double dis;
    double line[][] = new double[screenWidth][3];
    for (int y = 0; y < screenHeight; ++y){
      for (int x = 0; x < screenWidth; ++x){
        ray.x = firstray.x + x * scrnx.x - y * scrny.x;
        ray.y = firstray.y + x * scrnx.y - y * scrny.y;
        ray.z = firstray.z + x * scrnx.z - y * scrny.z;
        
        ray.normalize();
        
        dis = intersect(null,scene.eye,ray,color);
        
        if (dis > 0) /* ray intersected object  */
        {
          line[x][0] = color.r;
          line[x][1] = color.g;
          line[x][2] = color.b;
        } else /* use background color  */
        {
          line[x][0] = scene.background.r;
          line[x][1] = scene.background.g;
          line[x][2] = scene.background.b;
        }
      }
    }
  }
  
  private double intersect(Shape source,Point3d pos,Point3d ray,Color color){
    Shape objhit;
    double s, ss;
    Point3d hit, normal;
    Surface surf;
    double u, v;

    objhit = null;
    ss = Double.MAX_VALUE;
    /* check for intersection of ray with all objects */
    for (Shape objtry : scene.shapes) {
      if (objtry != source) /* don't try source */
      {
        s = objtry.intersect(pos, ray);
        /* keep track of closest intersection */
        if ((s > 0.0) && (s < ss)) {
          objhit = objtry;
          ss = s;
        }
      }
    }
    if (objhit == null)
      return 0.0; /* ray hit no objects */

    /* find point of intersection */
    
    hit = new Point3d(pos.x + ss*ray.x,
                      pos.y + ss*ray.y,
                      pos.z + ss*ray.z);
    /* find normal */
    normal = objhit.norm(hit);
    /* find color at point of intersection */
    shade(hit, ray, normal, objhit, color);
    surf = objhit.surface;
    /*
    if (surf->pattern == 1) {
      (*objpat[object[objhit].objtyp])(hit, &object[objhit], &u, &v);
      checkerbrd(u, v, 10, 10, color);
    }*/
    return (ss);
  }
  
  public void shade(Point3d pos, Point3d ray, Point3d nrm, Shape obj,
      Color color) {
    double k,bright, spec, diffuse;
    double nr, alpha, beta, cos_t;
    Surface surf = obj.surface;
    Point3d refl = new Point3d();
    Point3d trans = new Point3d();
    Color newcolor = new Color();
    double dist;

    /* ambient light contribution */
    color.r = surf.ambient.r;
    color.g = surf.ambient.g;
    color.b = surf.ambient.b;

    /* calculate reflected ray */
    k = -2.0 * Point3d.dotp(ray, nrm);
    refl.x = k * nrm.x + ray.x;
    refl.y = k * nrm.y + ray.y;
    refl.z = k * nrm.z + ray.z;
    if ((surf.reflection > 0.0) && (level <= depth)) {
      level++;
      dist = intersect(obj, pos, refl, newcolor);
      if (dist > 0) {
        color.r += (newcolor.r * surf.reflection);
        color.g += (newcolor.g * surf.reflection);
        color.b += (newcolor.b * surf.reflection);
      } else {
        color.r += (scene.background.r * surf.reflection);
        color.g += (scene.background.g * surf.reflection);
        color.b += (scene.background.b * surf.reflection);
      }
      level--;
    }

    if ((surf.transparency > 0) && (level <= depth)) {
      level++;
      nr = 1. / surf.medium;
      if (nr > 1.0)
        nr = 1.0;
      alpha = nr;
      cos_t = -Point3d.dotp(ray, nrm);
      beta = nr * cos_t - Math.sqrt(1. - nr * nr * (1. - cos_t * cos_t));
      trans.x = alpha * ray.x + beta * nrm.x;
      trans.y = alpha * ray.y + beta * nrm.y;
      trans.z = alpha * ray.z + beta * nrm.z;
      dist = intersect(obj, pos, trans, newcolor);
      if (dist > 0) {
        color.r = color.r * (1. - surf.transparency) + (newcolor.r * surf.transparency);
        color.g = color.g * (1. - surf.transparency) + (newcolor.g * surf.transparency);
        color.b = color.b * (1. - surf.transparency) + (newcolor.b * surf.transparency);
      } else {
        color.r = color.r * (1. - surf.transparency) + (scene.background.r * surf.transparency);
        color.g = color.g * (1. - surf.transparency) + (scene.background.g * surf.transparency);
        color.b = color.b * (1. - surf.transparency) + (scene.background.b * surf.transparency);
      }
      level--;
    }

    Point3d ltray = new Point3d();
    for (Light light : scene.lights) {
      /* get ray to light */
      lightray(light, pos, ltray);
      diffuse = Point3d.dotp(nrm, ltray);
      if (diffuse > 0.0) {
        /* object faces light, add diffuse */
        bright = brightness(obj, light, pos, ltray);
        diffuse *= bright;
        color.r += surf.diffuse.r * diffuse;
        color.g += surf.diffuse.g * diffuse;
        color.b += surf.diffuse.b * diffuse;

        spec = Point3d.dotp(refl, ltray);
        if (spec > 0) {
          /* highlight is here, add specular */
          spec = bright * Math.pow(spec, surf.coefficient);
          color.r += surf.specular.r * spec;
          color.g += surf.specular.g * spec;
          color.b += surf.specular.b * spec;
        }
      }
    }

  }


/* lightray() calculates vector from object point objpos to light source */
/* number lnum                                                           */
 
  private static void lightray(Light light, Point3d objpos, Point3d lray) {
    lray.x = light.position.x - objpos.x;
    lray.y = light.position.y - objpos.y;
    lray.z = light.position.z - objpos.z;
    lray.normalize();
  }
  
  /* brightness() returns intensity of light source # lnum as seen from  */
  /* position pos along ray direction -- Trivial for simple tracer, but  */
  /* is a shadow feeler for the recursive ray tracer. Will need to make  */
  /* calls to object intersection routines. If ray is blocked, returned  */
  /* brightness is zero.                                                 */

  private double brightness(Shape source, Light light, Point3d pos, Point3d ray) {
    int i;
    double t;
    double ss, s;
    t = (light.position.x - pos.x) / ray.x;
    if (t == 0)
      t = (light.position.y - pos.y) / ray.y;
    if (t == 0)
      t = (light.position.z - pos.z) / ray.z;

    ss = Double.MAX_VALUE;
    /* check for intersection of ray with all objects */
    for (Shape obj : scene.shapes) {
      if (obj != source) /* don't try source */
      {
        s = obj.intersect(pos, ray);
        /* keep track of closest intersection */
        if ((s > 0.0) && (s < ss))
          ss = s;
      }
    }
    if (ss < t)
      return 0.0;
    return (light.brightness);
  }
  
  public static void main(String[] args) throws Exception{
    Scene scene = Scene.makeScene();
    Renderer renderer = new Renderer(320,200,scene);
    renderer.initScene();
    renderer.render();
  }
}
