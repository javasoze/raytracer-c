package raytracerj;

import java.util.LinkedList;
import java.util.List;

public class Scene {
  public Point3d eye;
  public Point3d look;
  public Point3d up;
  public Color background;
  public List<Light> lights;
  public List<Shape> shapes;
  public double hfov;
  public double vfov;
  
  public Scene(){
    look = new Point3d(0,0,0);
    up = new Point3d(0,1,0);
    hfov = 60.0;
    vfov = 50.0;
    lights = new LinkedList<Light>();
    shapes = new LinkedList<Shape>();
    background = new Color(255,255,255);
  }
  
  public void addLight(Point3d light){
    this.addLight(light);
  }
  
  public void addShape(Shape shape){
    this.shapes.add(shape);
  }
  
  public static Scene makeScene(){
    Scene scene = new Scene();
    
    scene.eye = new Point3d(100,0,100);
    scene.look = new Point3d(0,0,0);
    
    Point3d light = new Point3d(100,150,50);
    scene.addLight(light);
    
    Surface redSurf = new Surface();
    redSurf.ambient = new Color(100,0,0);
    redSurf.diffuse = new Color(100,0,0);
    redSurf.specular = new Color(190,190,190);
    redSurf.coefficient = 30;
    redSurf.transparency = 0.8;
    redSurf.medium = 1.0;
    redSurf.pattern = 0.0;
    Sphere sph = new Sphere(redSurf,new Point3d(0,0,0),40);
    scene.addShape(sph);
    return scene;
  }
}
