package raytracerj;

public abstract class Shape {

  public final Surface surface;
  
  public Shape(Surface surface){
    this.surface = surface;
  }
  
  public abstract Point3d norm(Point3d pos);
  public abstract double intersect(Point3d pos, Point3d ray);
}
