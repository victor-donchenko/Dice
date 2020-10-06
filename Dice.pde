double interpolate(
  double before,
  double after,
  double progress
) {
  return before + progress * (after - before);
}

double rand_uniform(
  double bottom,
  double top
) {
  return bottom + Math.random() * (top - bottom);
}

int rand_int(
  int bottom,
  int top
) {
  return bottom + (int)(Math.random() * (top - bottom));
}

class Vector2d {
  public double x;
  public double y;

  Vector2d(double i_x, double i_y) {
    x = i_x;
    y = i_y;
  }
  
  Vector2d(Vector2d other) {
    this(other.x, other.y);
  }
}

class UniformRandomGenerator {
  public double bottom;
  public double top;
  
  UniformRandomGenerator(double i_bottom, double i_top) {
    bottom = i_bottom;
    top = i_top;
  }
  
  double get_random_number() {
    return rand_uniform(bottom, top);
  }
}

class Vector3d {
  public double x;
  public double y;
  public double z;

  Vector3d(double i_x, double i_y, double i_z) {
    x = i_x;
    y = i_y;
    z = i_z;
  }
  
  Vector3d(Vector3d other) {
    this(other.x, other.y, other.z);
  }
  
  void assign(Vector3d other) {
    x = other.x;
    y = other.y;
    z = other.z;
  }
  
  Vector3d plus(Vector3d other) {
    return new Vector3d(
      x + other.x,
      y + other.y,
      z + other.z
    );
  }
  
  Vector3d minus(Vector3d other) {
    return new Vector3d(
      x - other.x,
      y - other.y,
      z - other.z
    );
  }
  
  Vector3d multiply(double other) {
    return new Vector3d(
      x * other,
      y * other,
      z * other
    );
  }
  
  void update_by(Vector3d other) {
    x += other.x;
    y += other.y;
    z += other.z;
  }
  
  void apply_matrix(Matrix3d matrix) {
    assign(multiply_matrix(matrix, this));
  }
    
  String toString() {
    StringBuilder out = new StringBuilder();
    out.append('(');
    out.append(x);
    out.append(", ");
    out.append(y);
    out.append(", ");
    out.append(z);
    out.append(')');
    return out.toString();
  }
}

Vector3d interpolate(
  Vector3d before,
  Vector3d after,
  double progress
) {
  return new Vector3d(
    interpolate(before.x, after.x, progress),
    interpolate(before.y, after.y, progress),
    interpolate(before.z, after.z, progress)
  );
}

class Matrix3d {
  public double elements[][];
  
  Matrix3d(double i_elements[][]) {
    elements = new double[3][3];
    arrayCopy(i_elements, elements);
  }
}

Vector3d multiply_matrix(Matrix3d mat, Vector3d v) {
  return new Vector3d(
    mat.elements[0][0] * v.x + mat.elements[0][1] * v.y + mat.elements[0][2] * v.z,
    mat.elements[1][0] * v.x + mat.elements[1][1] * v.y + mat.elements[1][2] * v.z,
    mat.elements[2][0] * v.x + mat.elements[2][1] * v.y + mat.elements[2][2] * v.z
  );
}

Matrix3d get_rotation_matrix(Vector3d angles) {
  double x = angles.x;
  double y = angles.y;
  double z = angles.z;
  double cx = Math.cos(x);
  double sx = Math.sin(x);
  double cy = Math.cos(y);
  double sy = Math.sin(y);
  double cz = Math.cos(z);
  double sz = Math.sin(z);
  double elements[][] = {
    { cz * cy, cz * sy * sx - sz * cx, cz * sy * cx + sz * sx },
    { sz * cy, sz * sy * sx + cz * cx, sz * sy * cx - cz * sx },
    { -sy, cy * sx, cy * cx }
  };
  return new Matrix3d(elements);
}

class CubePosition {
  public Vector3d center;
  public Vector3d ang_displacement;
  
  CubePosition(
    Vector3d i_center,
    Vector3d i_ang_displacement
  ) {
    center = i_center;
    ang_displacement = i_ang_displacement;
  }
  
  CubePosition(CubePosition other) {
    this(new Vector3d(other.center), new Vector3d(other.ang_displacement));
  }
  
  String toString() {
    StringBuilder out = new StringBuilder();
    out.append("{center: ");
    out.append(center);
    out.append(", ang_displacement: ");
    out.append(ang_displacement);
    out.append('}');
    return out.toString();
  }
}

CubePosition interpolate(
  CubePosition before,
  CubePosition after,
  double progress
) {
  return new CubePosition(
    interpolate(
      before.center,
      after.center,
      progress
    ),
    interpolate(
      before.ang_displacement,
      after.ang_displacement,
      progress
    )
  );
}

class CubePath {
  public ArrayList<CubePosition> positions;
  public ArrayList<Double> times;
  
  CubePath() {
    positions = new ArrayList<CubePosition>();
    times = new ArrayList<Double>();
  }
  
  void add_waypoint(CubePosition pos, double time) {
    positions.add(pos);
    times.add(time);
  }
  
  void invert() {
    if (positions.size() == 0) {
      return;
    }
    double latest_time = times.get(times.size() - 1);
    for (int i = 0; i < positions.size() / 2; ++i) {
      int j = positions.size() - i - 1;
      CubePosition temp1 = positions.get(i);
      double temp2 = times.get(i);
      positions.set(i, positions.get(j));
      times.set(i, latest_time - times.get(j));
      positions.set(j, temp1);
      times.set(j, latest_time - temp2);
    }
  }
  
  String toString() {
    StringBuilder out = new StringBuilder();
    out.append("# positions: ");
    out.append(positions.size());
    out.append('\n');
    for (int i = 0; i < positions.size(); ++i) {
      out.append("  ");
      out.append(times.get(i));
      out.append(": ");
      out.append(positions.get(i));
      out.append('\n');
    }
    return out.toString();
  }
}

CubePath get_static_path(CubePosition static_position) {
  CubePath ret = new CubePath();
  ret.add_waypoint(static_position, 0);
  return ret;
}

class CubeMotionState {
  public int stage;
  public double time_until_next_stage; 
  public CubePosition position;
  public Vector3d velocity;
  public Vector3d ang_velocity;
  public double extra_data1;
  public double extra_data2;
  
  CubeMotionState(
    int i_stage,
    double i_time_until_next_stage,
    CubePosition i_position,
    Vector3d i_velocity,
    Vector3d i_ang_velocity,
    double i_extra_data1,
    double i_extra_data2
  ) {
    stage = i_stage;
    time_until_next_stage = i_time_until_next_stage;
    position = i_position;
    velocity = i_velocity;
    ang_velocity = i_ang_velocity;
    extra_data1 = i_extra_data1;
    extra_data2 = i_extra_data2;
  }
}

class Cube /* extends Shape */ {
  private double side_length;
  private CubePosition base_position;
  private CubePath path;
  private double current_time;
  private int last_path_index;

  Cube(double i_side_length, Vector3d i_base_center) {
    side_length = i_side_length;
    base_position = new CubePosition(
      i_base_center,
      new Vector3d(0, 0, 0)
    );
    path = get_static_path(base_position);
    current_time = 0;
    last_path_index = 0;
  }
  
  void reset() {
    path = get_static_path(base_position);
    current_time = 0;
    last_path_index = 0;
  }

  double get_side_length() {
    return side_length;
  }
  
  CubePosition get_base_position() {
    return base_position;
  }
  
  double get_current_time() {
    return current_time;
  }
  
  void add_time(double time_diff) {
    current_time += time_diff;
  }
  
  CubePath get_path() {
    return path;
  }
  
  void set_path(CubePath new_path) {
    path = new_path;
  }
  
  CubePosition get_current_position() {
    while (last_path_index < path.positions.size() - 1
            && current_time >= path.times.get(last_path_index + 1)) {
      ++last_path_index;
    }
    if (last_path_index >= path.positions.size() - 1) {
      return path.positions.get(path.positions.size() - 1);
    }
    CubePosition before = path.positions.get(last_path_index);
    CubePosition after = path.positions.get(last_path_index + 1);
    double before_time = path.times.get(last_path_index);
    double after_time = path.times.get(last_path_index + 1);
    double progress = (current_time - before_time) / (after_time - before_time);
    return interpolate(before, after, progress);
  }
  
  private void draw_face(int side_num) {
    final double circle_size = 1.0 / 5;
    final double circle_gap = 1.0 / 2;
    final Vector2d circle_centers[][] = {
      {},
      {new Vector2d(0, 0)},
      {new Vector2d(-1, -1), new Vector2d(1, 1)},
      {new Vector2d(-1, -1), new Vector2d(0, 0), new Vector2d(1, 1)},
      {new Vector2d(-1, -1), new Vector2d(-1, 1), new Vector2d(1, 1),
        new Vector2d(1, -1)},
      {new Vector2d(-1, -1), new Vector2d(-1, 1), new Vector2d(1, 1),
        new Vector2d(1, -1), new Vector2d(0, 0)},
      {new Vector2d(-1, -1), new Vector2d(-1, 0), new Vector2d(-1, 1),
        new Vector2d(1, -1), new Vector2d(1, 0), new Vector2d(1, 1)}
    };
    rect(-1, -1, 1, 1);
    for (int i = 0; i < side_num; ++i) {
      Vector2d circle_center = circle_centers[side_num][i];
      ellipse(
        (float)(circle_gap * circle_center.x),
        (float)(circle_gap * circle_center.y),
        (float)circle_size,
        (float)circle_size
      );
    }
  }
  
  void show() { 
    CubePosition position = get_current_position();
    
    pushMatrix();
    translate(
      (float)position.center.x,
      (float)position.center.y,
      (float)position.center.z
    );
    rotateZ((float)position.ang_displacement.z);
    rotateY((float)position.ang_displacement.y);
    rotateX((float)position.ang_displacement.x);
    scale((float)(side_length / 2));
    
    fill(color(0xff, 0xff, 0xff));
    stroke(color(0x00, 0x00, 0x00));
    strokeWeight(0.1);
    rectMode(CORNERS);
    ellipseMode(RADIUS);
    
    // close face (1)
    pushMatrix();
    translate(0, 0, 1);
    draw_face(1);
    popMatrix();
    
    // far face (6)
    pushMatrix();
    translate(0, 0, -1);
    draw_face(6);
    popMatrix();
    
    // left face (2)
    pushMatrix();
    translate(-1, 0, 0);
    rotateY(PI/2);
    draw_face(2);
    popMatrix();
    
    // right face (5)
    pushMatrix();
    translate(1, 0, 0);
    rotateY(-PI/2);
    draw_face(5);
    popMatrix();
    
    // top face
    pushMatrix();
    translate(0, 1, 0);
    rotateX(PI/2);
    draw_face(3);
    popMatrix();
    
    // bottom face
    pushMatrix();
    translate(0, -1, 0);
    rotateX(-PI/2);
    draw_face(4);
    popMatrix();
    
    popMatrix();
  }
  
  void roll() {
    base_position.ang_displacement.x = (PI / 2) * rand_int(0, 4);
    base_position.ang_displacement.y = (PI / 2) * rand_int(0, 4);
    base_position.ang_displacement.z = (PI / 2) * rand_int(0, 4);
  }
}

Cube cubes[] = new Cube[3 * 3];
static final double cube_side_length = 60;
static final double cube_gap = 150;
static final double margin = 90;
static final double dice_display_width = 540;

void init_cubes() {
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) {
      Vector3d base_position = new Vector3d(
        margin + cube_side_length / 2 + cube_gap * c,
        margin + cube_side_length / 2 + cube_gap * r,
        cube_side_length / 2
      );
      cubes[3 * r + c] = new Cube(
        cube_side_length,
        base_position
      );
    }
  }
}

final double BUFFER = 0.001;
final double gravity = 0.001;
final double friction = 0.005;
final double bounce_factor = 1.2;
final UniformRandomGenerator gen_slide_time
  = new UniformRandomGenerator(250, 500);
final UniformRandomGenerator gen_vel_argument
  = new UniformRandomGenerator(0, 2 * PI);
final UniformRandomGenerator gen_upward_velocity
  = new UniformRandomGenerator(0.2, 0.3);
final UniformRandomGenerator gen_ang_vel_comp
  = new UniformRandomGenerator(0, PI / 1000);

void plan_cube_paths() {
  //System.out.println("plan_cube_paths()");
  
  CubeMotionState motion_state_array[] = new CubeMotionState[cubes.length];
  CubePath path_array[] = new CubePath[cubes.length];
  
  for (int i = 0; i < cubes.length; ++i) {
    Cube cube = cubes[i];
    
    path_array[i] = new CubePath();
    path_array[i].add_waypoint(
      new CubePosition(cube.get_base_position()),
      0
    );
    
    motion_state_array[i] = new CubeMotionState(
      -1,
      0,
      new CubePosition(cube.get_base_position()),
      new Vector3d(0, 0, 0),
      new Vector3d(0, 0, 0),
      0,
      0
    );
  }
  
  //System.out.println("finished creating paths and motion states");
  
  final double time_delta = 10;
  double time_elapsed = 0;
  int num_finished = 0;
  
  while (num_finished < cubes.length) {
    for (int i = 0; i < cubes.length; ++i) {
      CubePosition position = new CubePosition(motion_state_array[i].position);
      path_array[i].add_waypoint(
        position,
        time_elapsed
      );
    }
    for (int i = 0; i < cubes.length; ++i) {
      motion_state_array[i].time_until_next_stage -= time_delta;
      if (motion_state_array[i].time_until_next_stage <= 0) {
        ++motion_state_array[i].stage;
        switch (motion_state_array[i].stage) {
          case 0:
            motion_state_array[i].time_until_next_stage = gen_slide_time.get_random_number();
            double vel_argument = gen_vel_argument.get_random_number();
            motion_state_array[i].velocity = new Vector3d(0, 0, 0);
            motion_state_array[i].ang_velocity = new Vector3d(0, 0, 0);
            motion_state_array[i].extra_data1 = 0;
            motion_state_array[i].extra_data2 = vel_argument;
            break;
          case 1:
            motion_state_array[i].time_until_next_stage = 3000;
            double upward_velocity = gen_upward_velocity.get_random_number();
            motion_state_array[i].velocity.z = upward_velocity;
            double ang_vel_x = gen_ang_vel_comp.get_random_number();
            double ang_vel_y = gen_ang_vel_comp.get_random_number();
            double ang_vel_z = gen_ang_vel_comp.get_random_number();
            motion_state_array[i].ang_velocity = new Vector3d(
              ang_vel_x,
              ang_vel_y,
              ang_vel_z
            );
            break;
          case 2:
            ++num_finished;
            break;
        }
      }
    }
    for (int i = 0; i < cubes.length; ++i) {
      motion_state_array[i].position.center.update_by(
        motion_state_array[i].velocity.multiply(
          time_delta
        )
      );
      motion_state_array[i].position.ang_displacement.update_by(
        motion_state_array[i].ang_velocity.multiply(
          time_delta
        )
      );
    }
    for (int i = 0; i < cubes.length; ++i) {
      CubePosition position = motion_state_array[i].position;
      double side_length = cubes[i].get_side_length();
      
      if (position.center.z + BUFFER <= side_length * (Math.sqrt(3) / 2)) {
        Vector3d vertices[] = new Vector3d[8];
        
        int j = 0;
        for (int x = -1; x <= 1; x += 2) {
          for (int y = -1; y <= 1; y += 2) {
            for (int z = -1; z <= 1; z += 2) {
              vertices[j] = new Vector3d(x, y, z);
              ++j;
            }
          }
        }
        
        double minimum_z = 2;
        for (int k = 0; k < 8; ++k) {
          vertices[k].apply_matrix(get_rotation_matrix(position.ang_displacement));
          if (vertices[k].z < minimum_z) {
            minimum_z = vertices[k].z;
          }
        }
        
        //System.out.println(minimum_z);
        
        minimum_z *= side_length / 2;
        
        //System.out.println(minimum_z);
        //System.out.println(position.center.z);
        
        if (position.center.z + minimum_z + BUFFER < 0) {
          position.center.z = -minimum_z;
          if (motion_state_array[i].velocity.z < 0) {
            motion_state_array[i].velocity.z *= -bounce_factor;
          }
        }
      }
    }
    for (int i = 0; i < cubes.length; ++i) {
      if (motion_state_array[i].stage == 0) {
        motion_state_array[i].extra_data1 += friction;
        double magnitude = motion_state_array[i].extra_data1;
        double argument = motion_state_array[i].extra_data2;
        motion_state_array[i].velocity.x = magnitude * Math.cos(argument);
        motion_state_array[i].velocity.y = magnitude * Math.sin(argument);
      }
      if (motion_state_array[i].stage == 1) {
        motion_state_array[i].velocity.z -= gravity * time_delta;
      }
    }
    time_elapsed += time_delta;
  }
  
  //System.out.println("finished calculating paths");
  
  for (int i = 0; i < cubes.length; ++i) {
    path_array[i].invert();
    cubes[i].set_path(path_array[i]);
  }
  
  //System.out.println("finished assigning paths");
}

void reset_cubes() {
  for (int i = 0; i < cubes.length; ++i) {
    cubes[i].reset();
  }
}

void roll_cubes() {
  for (int i = 0; i < cubes.length; ++i) {
    cubes[i].roll();
  }
}

void show_cubes() {
  for (int i = 0; i < cubes.length; ++i) {
    cubes[i].show();
  }
}

int last_draw_time;
boolean suspend_draw;
boolean in_draw;

void setup() {
  size(540, 540, P3D);
  init_cubes();
  last_draw_time = 0;
  suspend_draw = false;
  in_draw = false;
}

void draw() {
  if (suspend_draw) {
    return;
  }
  
  in_draw = true;
  background(color(0xff, 0xff, 0xff));
  show_cubes();
  
  int current_time = millis();
  for (int i = 0; i < cubes.length; ++i) {
    cubes[i].add_time(current_time - last_draw_time);
  }
  
  last_draw_time = current_time;
  in_draw = false;
}

void mousePressed() {
  suspend_draw = true;
  while (in_draw) {}
  roll_cubes();
  reset_cubes();
  plan_cube_paths();
  //System.out.println(cubes[0].get_path());
  suspend_draw = false;
}
