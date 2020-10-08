double interpolate_double(
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
  
  double magnitude_sq() {
    return x * x + y * y + z * z;
  }
  
  double magnitude() {
    return Math.sqrt(magnitude_sq());
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

Vector3d interpolate_Vector3d(
  Vector3d before,
  Vector3d after,
  double progress
) {
  return new Vector3d(
    interpolate_double(before.x, after.x, progress),
    interpolate_double(before.y, after.y, progress),
    interpolate_double(before.z, after.z, progress)
  );
}

class Matrix3d {
  public double elements[][];
  
  Matrix3d(double i_elements[][]) {
    elements = new double[3][3];
    arrayCopy(i_elements, elements);
  }
  
    Matrix3d(
    double a00, double a01, double a02,
    double a10, double a11, double a12,
    double a20, double a21, double a22
  ) {
    double[][] i_elements = {
      {a00, a01, a02},
      {a10, a11, a12},
      {a20, a21, a22}
    };
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

Matrix3d multiply_matrices(Matrix3d a, Matrix3d b) {
  Vector3d[] out_columns = new Vector3d[3];
  for (int c = 0; c < 3; ++c) {
    out_columns[c] = multiply_matrix(
      a,
      new Vector3d(
        b.elements[0][c],
        b.elements[1][c],
        b.elements[2][c]
      )
    );
  }
  
  double[][] elements = new double[3][3];
  for (int c = 0; c < 3; ++c) {
    elements[0][c] = out_columns[c].x;
    elements[1][c] = out_columns[c].y;
    elements[2][c] = out_columns[c].z;
  }
  
  return new Matrix3d(elements);
}

Matrix3d get_x_rotation_matrix(double theta) {
  double ct = Math.cos(theta);
  double st = Math.sin(theta);
  return new Matrix3d(
    1,  0,   0,
    0, ct, -st,
    0, st,  ct
  );
}

Matrix3d get_y_rotation_matrix(double theta) {
  double ct = Math.cos(theta);
  double st = Math.sin(theta);
  return new Matrix3d(
     ct, 0, st,
      0, 1,  0,
    -st, 0, ct
  );
}

Matrix3d get_z_rotation_matrix(double theta) {
  double ct = Math.cos(theta);
  double st = Math.sin(theta);
  return new Matrix3d(
     ct, -st, 0,
     st,  ct, 0,
      0,   0, 1
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

Matrix3d get_inverted_rotation_matrix(Vector3d angles) {
  return multiply_matrices(
    get_x_rotation_matrix(-angles.x),
    multiply_matrices(
      get_y_rotation_matrix(-angles.y),
      get_z_rotation_matrix(-angles.z)
    )
  );
}

class AffineTransform {
  public Matrix3d linear;
  public Vector3d translate;

  AffineTransform(Matrix3d i_linear, Vector3d i_translate) {
    linear = i_linear;
    translate = i_translate;
  }
}

Vector3d do_affine(AffineTransform affine, Vector3d v) {
  return multiply_matrix(affine.linear, v).plus(affine.translate);
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

CubePosition interpolate_CubePosition(
  CubePosition before,
  CubePosition after,
  double progress
) {
  return new CubePosition(
    interpolate_Vector3d(
      before.center,
      after.center,
      progress
    ),
    interpolate_Vector3d(
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
  
  double get_total_time() {
    return times.get(times.size() - 1);
  }
  
  void invert() {
    if (positions.size() == 0) {
      return;
    }
    double latest_time = get_total_time();
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

Vector3d[] get_unit_circle_lines() {
  final int num_vertices = 20;

  Vector3d[] vertices = new Vector3d[num_vertices];
  for (int i = 0; i < num_vertices; ++i) {
    double angle = 2*PI * i/num_vertices;
    vertices[i] = new Vector3d(Math.cos(angle), Math.sin(angle), 0);
  }

  Vector3d[] lines_vertices = new Vector3d[num_vertices * 2];
  for (int i = 0; i < num_vertices; ++i) {
    lines_vertices[2*i] = vertices[i];
    lines_vertices[2*i+1] = vertices[(i+1) % num_vertices];
  }
  
  return lines_vertices;
}

Vector3d[] get_dice_quads() {
  Vector3d vertices[] = {
    new Vector3d(-1, -1, -1), // 0
    new Vector3d(-1, -1,  1), // 1
    new Vector3d(-1,  1, -1), // 2
    new Vector3d(-1,  1,  1), // 3
    new Vector3d( 1, -1, -1), // 4
    new Vector3d( 1, -1,  1), // 5
    new Vector3d( 1,  1, -1), // 6
    new Vector3d( 1,  1,  1)  // 7
  };

  int faces[][] = {
    {0, 1, 3, 2},
    {0, 1, 5, 4},
    {0, 2, 6, 4},
    {7, 6, 4, 5},
    {7, 6, 2, 3},
    {7, 5, 1, 3}
  };

  Vector3d[] out = new Vector3d[4 * faces.length];

  for (int i = 0; i < faces.length; ++i) {
    for (int j = 0; j < 4; ++j) {
      out[4*i + j] = vertices[faces[i][j]];
    }
  }

  return out;
}

Vector3d[] get_dice_lines() {
  final Vector3d[] unit_circle_lines = get_unit_circle_lines();

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

  final AffineTransform affine_transforms[] = {
    new AffineTransform(new Matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0), new Vector3d(0, 0, 0)), // 0
    new AffineTransform(get_rotation_matrix(new Vector3d(0, 0, 0)), new Vector3d(0, 0, 1)), // 1
    new AffineTransform(get_rotation_matrix(new Vector3d(0, PI/2, 0)), new Vector3d(-1, 0, 0)), // 2
    new AffineTransform(get_rotation_matrix(new Vector3d(PI/2, 0, 0)), new Vector3d(0, 1, 0)), // 3
    new AffineTransform(get_rotation_matrix(new Vector3d(-PI/2, 0, 0)), new Vector3d(0, -1, 0)), // 4
    new AffineTransform(get_rotation_matrix(new Vector3d(0, -PI/2, 0)), new Vector3d(1, 0, 0)), // 5
    new AffineTransform(get_rotation_matrix(new Vector3d(0, 0, 0)), new Vector3d(0, 0, -1)), // 6
  };

  ArrayList<Vector3d> lines = new ArrayList<Vector3d>();

  for (int s = 1; s <= 6; ++s) {
    Vector2d[] centers = circle_centers[s];
    AffineTransform affine = affine_transforms[s];
    for (int ci = 0; ci < centers.length; ++ci) {
      Vector2d center = centers[ci];
      Vector3d center_diff = new Vector3d(
        center.x * circle_gap,
        center.y * circle_gap,
        0
      );
      for (int vi = 0; vi < unit_circle_lines.length; ++vi) {
        Vector3d vertex = unit_circle_lines[vi];
        vertex = vertex.multiply(circle_size);
        vertex = vertex.plus(center_diff);
        vertex = do_affine(affine, vertex);
        lines.add(vertex);
      }
    }
  }

  Vector3d[] lines_array = new Vector3d[lines.size()];
  for (int i = 0; i < lines.size(); ++i) {
    lines_array[i] = lines.get(i);
  }
  
  return lines_array;
}

final Vector3d[] dice_quads = get_dice_quads();
final Vector3d[] dice_lines = get_dice_lines();

class Cube {
  private double side_length;
  private CubePosition base_position;
  private CubePath path;
  private double current_time;
  private int last_path_index;
  private int roll_num;

  Cube(double i_side_length, Vector3d i_base_center) {
    side_length = i_side_length;
    base_position = new CubePosition(
      i_base_center,
      new Vector3d(0, 0, 0)
    );
    path = get_static_path(base_position);
    current_time = 0;
    last_path_index = 0;
    roll_num = 1;
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
  
  int get_roll_num() {
    return roll_num;
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
    return interpolate_CubePosition(before, after, progress);
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

    beginShape(QUADS);
    for (int i = 0; i < dice_quads.length; ++i) {
      Vector3d v = dice_quads[i];
      vertex((float)v.x, (float)v.y, (float)v.z);
    }
    endShape();

    beginShape(LINES);
    for (int i = 0; i < dice_lines.length; ++i) {
      Vector3d v = dice_lines[i];
      vertex((float)v.x, (float)v.y, (float)v.z);
    }
    endShape();
    
    popMatrix();
  }
  
  void roll() {
    roll_num = rand_int(1, 7);
    int x_displ_mult = 0;
    int y_displ_mult = 0;
    switch (roll_num) {
      case 1:
        x_displ_mult = 0;
        y_displ_mult = 0;
        break;
      case 2:
        x_displ_mult = 0;
        y_displ_mult = 1;
        break;
      case 3:
        x_displ_mult = 1;
        y_displ_mult = 0;
        break;
      case 4:
        x_displ_mult = -1;
        y_displ_mult = 0;
        break;
      case 5:
        x_displ_mult = 0;
        y_displ_mult = -1;
        break;
      case 6:
        x_displ_mult = 2;
        y_displ_mult = 0;
        break;
    }
    base_position.ang_displacement.x = (PI / 2) * x_displ_mult;
    base_position.ang_displacement.y = (PI / 2) * y_displ_mult;
    base_position.ang_displacement.z = rand_uniform(0, 2 * PI);
  }
}

class Planner {
  private final double BUFFER = 0.001;
  private final double gravity = 0.001;
  private final double friction = 0.005;
  private final double bounce_factor = 1.2;
  private final UniformRandomGenerator gen_slide_time
    = new UniformRandomGenerator(250, 500);
  private final UniformRandomGenerator gen_vel_argument
    = new UniformRandomGenerator(0, 2 * PI);
  private final UniformRandomGenerator gen_upward_velocity
    = new UniformRandomGenerator(0.2, 0.3);
  private final UniformRandomGenerator gen_ang_vel_comp
    = new UniformRandomGenerator(0, PI / 1000);
  private final double time_delta = 10;

  private Cube cube_array[];
  private CubeMotionState motion_state_array[];
  private CubePath path_array[];
  private double time_elapsed;
  private int num_finished;
  private HashMap<Integer, Integer> prev_collisions;
  
  Planner(Cube[] i_cube_array) {
    cube_array = i_cube_array;
    motion_state_array = null;
    path_array = null;
    time_elapsed = 0;
    num_finished = 0;
    prev_collisions = null;
  }
  
  void initialize() {
    motion_state_array = new CubeMotionState[cube_array.length];
    path_array = new CubePath[cube_array.length];
    
    for (int i = 0; i < cube_array.length; ++i) {
      Cube cube = cube_array[i];
      
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
  }
  
  void record_waypoints() {
    for (int i = 0; i < cube_array.length; ++i) {
      CubePosition position = new CubePosition(motion_state_array[i].position);
      path_array[i].add_waypoint(
        position,
        time_elapsed
      );
    }
  }
  
  void manage_stages() {
    for (int i = 0; i < cube_array.length; ++i) {
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
  }
  
  void move_positions() {
    for (int i = 0; i < cube_array.length; ++i) {
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
  }

  boolean precheck_collision_with_ground(int i) {
    CubePosition position = motion_state_array[i].position;
    double side_length = cube_array[i].get_side_length();
    return position.center.z + BUFFER <= side_length * (Math.sqrt(3) / 2);
  }
  
  boolean check_and_handle_collision_with_ground(int i) {
    CubePosition position = motion_state_array[i].position;
    double side_length = cube_array[i].get_side_length();
    
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
          
    minimum_z *= side_length / 2;
            
    if (position.center.z + minimum_z + BUFFER < 0) {
      position.center.z = -minimum_z;
      if (motion_state_array[i].velocity.z < 0) {
        motion_state_array[i].velocity.z *= -bounce_factor;
      }
      return true;
    }
    else {
      return false;
    }
  }

  boolean precheck_collision_of_cubes(int i1, int i2) {
    Vector3d center1 = motion_state_array[i1].position.center;
    Vector3d center2 = motion_state_array[i2].position.center;
    double side_length1 = cube_array[i1].side_length;
    double side_length2 = cube_array[i2].side_length;
    return center1.minus(center2).magnitude() < (side_length1 + side_length2) * Math.sqrt(3) / 2; 
  }

  double[] elastic_results(double m1, double v1, double m2, double v2) {
    double m_sum = m1 + m2;
    double m_diff = m1 - m2;
    double w1 = (m_diff * v1 + 2 * m2 * v2) / m_sum;
    double w2 = (2 * m1 * v1 - m_diff * v2) / m_sum;
    double[] out = {w1, w2};
    return out;
  }
  
  void elastic_rebound(int i1, int i2) {
    double mass1 = Math.pow(cube_array[i1].side_length, 3);
    Vector3d velocity1 = motion_state_array[i1].velocity;
    double mass2 = Math.pow(cube_array[i2].side_length, 3);
    Vector3d velocity2 = motion_state_array[i2].velocity;
    
    double x_results[] = elastic_results(mass1, velocity1.x, mass2, velocity2.x);
    double y_results[] = elastic_results(mass1, velocity1.y, mass2, velocity2.y);
    double z_results[] = elastic_results(mass1, velocity1.z, mass2, velocity2.z);
    
    motion_state_array[i1].velocity = new Vector3d(
      x_results[0],
      y_results[0],
      z_results[0]
    );
    
    motion_state_array[i2].velocity = new Vector3d(
      x_results[1],
      y_results[1],
      z_results[1]
    );
  }

  boolean check_and_handle_collision_to_cube(int i1, int i2) {
    double side_length1 = cube_array[i1].side_length;
    double side_length2 = cube_array[i2].side_length;
    CubePosition position1 = motion_state_array[i1].position;
    CubePosition position2 = motion_state_array[i2].position;
    CubePosition position_diff = new CubePosition(
      position2.center.minus(position1.center),
      position2.ang_displacement.minus(position1.ang_displacement)
    );
    
    Vector3d vertices[] = new Vector3d[64];
      
    int j = 0;
    for (double x = -1; x <= 1 + BUFFER; x += 2.0/3) {
      for (double y = -1; y <= 1 + BUFFER; y += 2.0/3) {
        for (double z = -1; z <= 1 + BUFFER; z += 2.0/3) {
          vertices[j] = (new Vector3d(x, y, z)).multiply(side_length2 / 2);
          ++j;
        }
      }
    }

    Matrix3d inverse_matrix = get_inverted_rotation_matrix(position1.ang_displacement);
    
    AffineTransform affine = new AffineTransform(
      multiply_matrices(
        inverse_matrix,
        get_rotation_matrix(position2.ang_displacement)
      ),
      multiply_matrix(
        inverse_matrix,
        position_diff.center
      )
    );
    
    for (j = 0; j < vertices.length; ++j) {
      vertices[j] = do_affine(affine, vertices[j]);
    }
    
    boolean has_been_collision = false;
    for (j = 0; j < vertices.length; ++j) {
      Vector3d v = vertices[j];
      double cube2_offset = side_length1 / 2;
      if (Math.abs(v.x) + BUFFER < cube2_offset
          && Math.abs(v.y) + BUFFER < cube2_offset
          && Math.abs(v.z) + BUFFER < cube2_offset) {
        has_been_collision = true;
        break;
      }
    }
    
    if (has_been_collision) {
      if (!prev_collisions.containsKey(cube_array.length * i1 + i2)) {
        elastic_rebound(i1, i2);
      }
      return true;
    }
    else {
      return false;
    }
  }
  
  void check_and_handle_collisions() {
    HashMap<Integer, Integer> current_collisions = new HashMap<Integer, Integer>();
    
    for (int i = 0; i < cube_array.length; ++i) {
      if (precheck_collision_with_ground(i)) {
        check_and_handle_collision_with_ground(i);
      }

      int i1 = i;
      for (int i2 = i + 1; i2 < cube_array.length; ++i2) {
        if (precheck_collision_of_cubes(i1, i2)) {
          boolean collision_happened = false;
          if (check_and_handle_collision_to_cube(i1, i2)) {
            collision_happened = true;
          }
          else if (check_and_handle_collision_to_cube(i2, i1)) {
            collision_happened = true;
          }
          if (collision_happened) {
            current_collisions.put(cube_array.length * i1 + i2, 1);
          }
        }
      }
    }
    
    prev_collisions = current_collisions;
  }
  
  void do_acceleration() {
    for (int i = 0; i < cube_array.length; ++i) {
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
  }
  
  double plan_cube_paths() {
    initialize();
    
    time_elapsed = 0;
    num_finished = 0;
    prev_collisions = new HashMap<Integer, Integer>();
    
    while (num_finished < cube_array.length) {
      record_waypoints();
      manage_stages();
      move_positions();
      check_and_handle_collisions();
      do_acceleration();
      
      time_elapsed += time_delta;
    }
    
    double longest_time = 0;
    
    for (int i = 0; i < cube_array.length; ++i) {
      if (path_array[i].get_total_time() > longest_time) {
        longest_time = path_array[i].get_total_time();
      }
      
      path_array[i].invert();
      cube_array[i].set_path(path_array[i]);
    }
    
    return longest_time;
  }
}

Cube cubes[] = new Cube[3*3];
Planner planner = new Planner(cubes);
String bottom_label_text;

final double cube_side_length = 60;
final double cube_gap = 150;
final double margin = 90;
final double dice_display_width = 540;
final double bottom_label_height = 60;
final double total_height = 600;

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

void init_bottom_label() {
  bottom_label_text = "Sum of all rolls: 9";
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

void change_bottom_label_text() {
  int roll_sum = 0;
  for (int i = 0; i < cubes.length; ++i) {
    roll_sum += cubes[i].get_roll_num();
  }
  bottom_label_text = "Sum of all rolls: " + roll_sum;
}

void show_cubes() {
  for (int i = 0; i < cubes.length; ++i) {
    cubes[i].show();
  }
}

void show_bottom_label() {
  rectMode(CORNER);
  fill(color(0x00, 0x00, 0x00));
  textSize(12);
  textAlign(CENTER);
  text(
    bottom_label_text,
    0, (int)dice_display_width,
    (float)dice_display_width, (float)bottom_label_height
  );
}

int last_draw_time;
boolean suspend_draw;
boolean in_draw;
double time_until_label_show;
int view;
final int num_views = 3;

void setup() {
  size(540, 600, P3D);
  background(color(0xff, 0xff, 0xff));
  
  init_cubes();
  init_bottom_label();
  last_draw_time = 0;
  suspend_draw = false;
  in_draw = false;
  time_until_label_show = 0;
  view = 1;
}

void view_change() {
  switch (view) {
    case 0:
      break;
    case 1:
      translate(0, 0, -280);
      rotateX(PI / 6);
      
      fill(color(0x80, 0x80, 0x80));
      rectMode(CORNERS);
      rect(
        0,
        0,
        540,
        540
      );
      break;
    case 2:
      stroke(color(0x80, 0x80, 0x80));
      strokeWeight(2);
      line(
        0,
        300,
        540,
        300
      );
      
      translate(0, 300, -300);
      rotateX(PI / 2);
      break;
  }
}

void draw() {
  if (suspend_draw) {
    return;
  }
  
  in_draw = true;
  
  background(color(0xff, 0xff, 0xff));
  
  pushMatrix();
  view_change();
  show_cubes();
  popMatrix();
  
  if (time_until_label_show <= 0) {
    show_bottom_label();
  }
  
  int current_time = millis();
  for (int i = 0; i < cubes.length; ++i) {
    cubes[i].add_time(current_time - last_draw_time);
  }
  if (time_until_label_show > 0) {
    time_until_label_show -= current_time - last_draw_time;
  }
  
  last_draw_time = current_time;
  in_draw = false;
}

void mousePressed() {
  suspend_draw = true;
  while (in_draw) {}
  
  if (mouseY > dice_display_width) {
    ++view;
    view %= num_views;
  }
  else {
    roll_cubes();
    change_bottom_label_text();
    reset_cubes();
    time_until_label_show = planner.plan_cube_paths();
  }
  
  suspend_draw = false;
}
