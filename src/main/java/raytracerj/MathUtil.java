package raytracerj;

public class MathUtil {

  public static final double radian(double degree){
    return degree * Math.PI / 180.0;
  }
  
  /* find roots of a quadratic polynomial */
  public static final boolean solveQuad(double a, double b, double c, double[] roots) {
    double discrim = b*b - 4. * a * c;

    if (a == 0)
      return false;
    if (discrim < 0.)
      return false;

    roots[0] = (-b + Math.sqrt(discrim)) / (2. * a);
    roots[1] = (-b - Math.sqrt(discrim)) / (2. * a);
    return true;
  }
  
  /**
   * Check to see if b is in [a,c]
   */
  public static final boolean between(double b, double a, double c) {
    if ((b <= c) && (b >= a))
      return true;
    return false;
  }
  
  public static double min_pos(double n1, double n2) {
    if ((n1 <= 0) && (n2 <= 0))
      return 0;
    if (n1 <= 0)
      return n2;
    if (n2 <= 0)
      return n1;
    if (n1 < n2)
      return n1;
    return n2;
  }
}
