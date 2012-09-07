package raytracerj;

public class Renderer {

  public static final int DEFAULT_DEPTH = 5;
  
  private int screenWidth;
  
  private int screenHeight;
  private int depth = DEFAULT_DEPTH;
  
  public Renderer(int screenWidth,int screenHeight){
    this.screenWidth = screenWidth;
    this.screenHeight = screenHeight;
  }
  
  public void initScene(){
    
  }
  
  public void render(Scene scene){
    Point3d scrnx = new Point3d();
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
        
        ray.noramlize();
        
        dis = intersect();
        
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
  
  private double intersect(){
    return 0.0;
  }
  
  public static void main(String[] args) {
    
  }
}
