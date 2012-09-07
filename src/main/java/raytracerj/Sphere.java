package raytracerj;

public class Sphere extends Shape {

  private Point3d center;
  private double radius;
  
  public Sphere(Surface surface,Point3d center,double radius){
    super(surface);
    this.center = center;
    this.radius = radius;
  }

  @Override
  public Point3d norm(Point3d pos) {
    Point3d norm = new Point3d(pos.x-center.x,pos.y-center.y,pos.z-center.z);
    norm.normalize();
    return norm;
  }

  @Override
  public double intersect(Point3d pos, Point3d ray) {
    double C1,C2,C3;
    double A,B;
    double C;
    C1 = radius*radius - center.x*center.x - center.y*center.y
        - center.z*center.z;
    C2 = C1 - pos.x*pos.x - pos.y*pos.y - pos.z*pos.z;

    C3 = C2 + 2. * pos.x * center.x + 2. * pos.y * center.y
        + 2. * pos.z * center.z;

    A = ray.x + ray.y*ray.y + ray.z*ray.z;
    B = 2.
        * (pos.x * ray.x - ray.x * center.x + pos.y * ray.y
            - ray.y * center.y + pos.z * ray.z - ray.z * center.z);
    C = -C3;
    double roots[] = new double[2];
    if (!MathUtil.solveQuad(A, B, C, roots))
      return Double.MAX_VALUE;
    return MathUtil.min_pos(roots[0], roots[1]);
  }
}
