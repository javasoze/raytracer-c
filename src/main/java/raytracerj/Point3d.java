package raytracerj;

public class Point3d {
  public double x=0,y=0,z=0;
  
  public final void noramlize(){
    double d = magnitude();
    if (d>0){
      x /= d;
      y /= d;
      z /= d;
    }
  }
  
  public final double magnitude(){
    return Math.sqrt(x * x + y * y + z * z);
  }
  
  public final void scale(double c){
    x *= c;
    y *= c;
    z *= c;
  }
  
  public final static double dotp(Point3d p1, Point3d p2){
    return (p1.x * p2.x + p1.y * p2.y + p1.z * p2.z);
  }
  
  public final static void crossp(Point3d a, Point3d b,Point3d result){
    result.x = (a.y * b.z) - (a.z * b.y);
    result.y = (a.z * b.x) - (a.x * b.z);
    result.z = (a.x * b.y) - (a.y * b.x);
  }
  
  public final static void copy(Point3d from,Point3d to){
    to.x = from.x;
    to.y = from.y;
    to.z = from.z;
  }
}

