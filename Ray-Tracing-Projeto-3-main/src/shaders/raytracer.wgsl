const THREAD_COUNT = 16;
const RAY_TMIN = 0.0001;
const RAY_TMAX = 100.0;
const PI = 3.1415927f;
const FRAC_1_PI = 0.31830987f;
const FRAC_2_PI = 1.5707964f;

@group(0) @binding(0)  
  var<storage, read_write> fb : array<vec4f>;

@group(0) @binding(1)
  var<storage, read_write> rtfb : array<vec4f>;

@group(1) @binding(0)
  var<storage, read_write> uniforms : array<f32>;

@group(2) @binding(0)
  var<storage, read_write> spheresb : array<sphere>;

@group(2) @binding(1)
  var<storage, read_write> quadsb : array<quad>;

@group(2) @binding(2)
  var<storage, read_write> boxesb : array<box>;

@group(2) @binding(3)
  var<storage, read_write> trianglesb : array<triangle>;

@group(2) @binding(4)
  var<storage, read_write> meshb : array<mesh>;

struct ray {
  origin : vec3f,
  direction : vec3f,
};

struct sphere {
  transform : vec4f,
  color : vec4f,
  material : vec4f,
};

struct quad {
  Q : vec4f,
  u : vec4f,
  v : vec4f,
  color : vec4f,
  material : vec4f,
};

struct box {
  center : vec4f,
  radius : vec4f,
  rotation: vec4f,
  color : vec4f,
  material : vec4f,
};

struct triangle {
  v0 : vec4f,
  v1 : vec4f,
  v2 : vec4f,
};

struct mesh {
  transform : vec4f,
  scale : vec4f,
  rotation : vec4f,
  color : vec4f,
  material : vec4f,
  min : vec4f,
  max : vec4f,
  show_bb : f32,
  start : f32,
  end : f32,
};

struct material_behaviour {
  scatter : bool,
  direction : vec3f,
};

struct camera {
  origin : vec3f,
  lower_left_corner : vec3f,
  horizontal : vec3f,
  vertical : vec3f,
  u : vec3f,
  v : vec3f,
  w : vec3f,
  lens_radius : f32,
};

struct hit_record {
  t : f32,
  p : vec3f,
  normal : vec3f,
  object_color : vec4f,
  object_material : vec4f,
  frontface : bool,
  hit_anything : bool,
};

fn ray_at(r: ray, t: f32) -> vec3f
{
  return r.origin + t * r.direction;
}

fn get_ray(cam: camera, uv: vec2f, rng_state: ptr<function, u32>) -> ray
{
  var rd = cam.lens_radius * rng_next_vec3_in_unit_disk(rng_state);
  var offset = cam.u * rd.x + cam.v * rd.y;
  return ray(cam.origin + offset, normalize(cam.lower_left_corner + uv.x * cam.horizontal + uv.y * cam.vertical - cam.origin - offset));
}

fn get_camera(lookfrom: vec3f, lookat: vec3f, vup: vec3f, vfov: f32, aspect_ratio: f32, aperture: f32, focus_dist: f32) -> camera
{
  var camera = camera();
  camera.lens_radius = aperture / 2.0;

  var theta = degrees_to_radians(vfov);
  var h = tan(theta / 2.0);
  var w = aspect_ratio * h;

  camera.origin = lookfrom;
  camera.w = normalize(lookfrom - lookat);
  camera.u = normalize(cross(vup, camera.w));
  camera.v = cross(camera.u, camera.w);

  camera.lower_left_corner = camera.origin - w * focus_dist * camera.u - h * focus_dist * camera.v - focus_dist * camera.w;
  camera.horizontal = 2.0 * w * focus_dist * camera.u;
  camera.vertical = 2.0 * h * focus_dist * camera.v;

  return camera;
}

fn envoriment_color(direction: vec3f, color1: vec3f, color2: vec3f) -> vec3f
{
  var unit_direction = normalize(direction);
  var t = 0.5 * (unit_direction.y + 1.0);
  var col = (1.0 - t) * color1 + t * color2;

  var sun_direction = normalize(vec3(uniforms[13], uniforms[14], uniforms[15]));
  var sun_color = int_to_rgb(i32(uniforms[17]));
  var sun_intensity = uniforms[16];
  var sun_size = uniforms[18];

  var sun = clamp(dot(sun_direction, unit_direction), 0.0, 1.0);
  col += sun_color * max(0, (pow(sun, sun_size) * sun_intensity));

  return col;
}

fn check_ray_collision(r: ray, max: f32) -> hit_record
{
    var spheresCount = i32(uniforms[19]);
    var quadsCount = i32(uniforms[20]);
    var boxesCount = i32(uniforms[21]);
    var trianglesCount = i32(uniforms[22]);
    var meshCount = i32(uniforms[27]);

    var closest = hit_record(max, vec3f(0.0), vec3f(0.0), vec4f(0.0), vec4f(0.0), false, false);

    // Temporary hit record for each object
    var temp_rec = hit_record(0.0, vec3f(0.0), vec3f(0.0), vec4f(0.0), vec4f(0.0), false, false);

    // Check spheres
    for (var i = 0; i < spheresCount; i = i + 1)
    {
        var sphere = spheresb[i];
        hit_sphere(sphere.transform.xyz, sphere.transform.w, r, &temp_rec, closest.t);

        if (temp_rec.hit_anything && temp_rec.t < closest.t)
        {
            closest = temp_rec;
            closest.object_color = sphere.color;
            closest.object_material = sphere.material;

            if (closest.object_material.x == 2.0) {
                var debug_detected_dialetic = true;
            }
        }
    }

    // Check quads
    for (var i = 0; i < quadsCount; i = i + 1)
    {
        var quad = quadsb[i];
        hit_quad(r, quad.Q, quad.u, quad.v, &temp_rec, closest.t);

        if (temp_rec.hit_anything && temp_rec.t < closest.t)
        {
            closest = temp_rec;
            closest.object_color = quad.color;
            closest.object_material = quad.material;
        }
    }

    // Check boxes
    for (var i = 0; i < boxesCount; i = i + 1)
    {
        var box_obj = boxesb[i];
        hit_box(r, box_obj.center.xyz, box_obj.radius.xyz, &temp_rec, closest.t);

        if (temp_rec.hit_anything && temp_rec.t < closest.t)
        {
            closest = temp_rec;
            closest.object_color = box_obj.color;
            closest.object_material = box_obj.material;
        }
    }

    // Check triangles
    for (var i = 0; i < trianglesCount; i = i + 1)
    {
        var tri = trianglesb[i];
        hit_triangle(r, tri.v0.xyz, tri.v1.xyz, tri.v2.xyz, &temp_rec, closest.t);

        if (temp_rec.hit_anything && temp_rec.t < closest.t)
        {
            closest = temp_rec;
            // Defina valores padrão se triângulos não têm materiais
            closest.object_color = vec4f(1.0);     // Cor padrão
            closest.object_material = vec4f(0.0, 0.0, 0.0, 1.0); // Default Lambertiano
        }
    }

    // Adicione lógica para meshes, se necessário.

    return closest;
}

fn lambertian(normal : vec3f, absorption: f32, random_sphere: vec3f, rng_state: ptr<function, u32>) -> material_behaviour
{
    // Calculate the scatter direction for Lambertian reflection
    var scatter_direction = normalize(normal + random_sphere);

    // Handle the case where scatter_direction is zero
    var epsilon = 1e-8;
    if (length(scatter_direction) < epsilon)
    {
        scatter_direction = normal;
    }

    // Return the material behaviour indicating that the ray should scatter
    return material_behaviour(true, scatter_direction);
}

fn metal(normal : vec3f, direction: vec3f, fuzz: f32, random_sphere: vec3f) -> material_behaviour
{
    // Reflect the incoming direction around the normal
    var reflected = reflect(normalize(direction), normal);

    // Add fuzziness to the reflection
    var scatter_direction = reflected + fuzz * random_sphere;

    // Normalize the scatter direction
    scatter_direction = normalize(scatter_direction);

    // Determine if the scattered ray is in the same hemisphere as the normal
    var scatter = dot(scatter_direction, normal) > 0.0;

    // Return the material behaviour indicating whether to scatter and the new direction
    return material_behaviour(scatter, scatter_direction);}

fn dielectric(normal : vec3f, r_direction: vec3f, refraction_index: f32, frontface: bool, random_sphere: vec3f, fuzz: f32, rng_state: ptr<function, u32>) -> material_behaviour
{  
    // Normalizando vetores
    var normalized_r_dir = normalize(r_direction);
    var normalized_normal = normalize(normal);

    // Índice de refração relativo
    var eta: f32;
    if (frontface) {
        eta = 1.0 / refraction_index;
    } else {
        eta = refraction_index;
    }

    // Cálculo de cos(theta) e sin(theta)
    var cos_theta = clamp(dot(-normalized_r_dir, normalized_normal), 0.0, 1.0);
    var sin_theta = sqrt(1.0 - cos_theta * cos_theta);

    // Verificação de reflexão total interna
    var cannot_refract = eta * sin_theta > 1.0;

    // Aproximação de Schlick para refvarância
    var r0 = (1.0 - eta) / (1.0 + eta);
    r0 = r0 * r0;
    var reflectance = r0 + (1.0 - r0) * pow(1.0 - cos_theta, 5.0);

    // Determinação da direção (reflexão ou refração)
    var direction: vec3f;
    if (cannot_refract || reflectance > rng_next_float(rng_state)) {
        direction = reflect(normalized_r_dir, normalized_normal); // Reflexão
    } else {
        direction = refract(normalized_r_dir, normalized_normal, eta); // Refração
    }

    // Retornando o comportamento do material
    return material_behaviour(true, normalize(direction));
}

fn emmisive(color: vec3f, light: f32) -> material_behaviour
{
    return material_behaviour(false, vec3f(0.0));
}

fn trace(r: ray, rng_state: ptr<function, u32>) -> vec3f
{
    var max_bounces = i32(uniforms[2]); // Maximum number of bounces
    var light = vec3f(0.0);             // Accumulated light
    var color = vec3f(1.0);             // Accumulated color
    var ray_current = r;                // Current ray being traced

    // Background colors from uniforms
    var background_color1 = int_to_rgb(i32(uniforms[11]));
    var background_color2 = int_to_rgb(i32(uniforms[12]));

    for (var bounce = 0; bounce < max_bounces; bounce = bounce + 1)
    {
        // Check for ray-object intersections
        var hit_rec = check_ray_collision(ray_current, RAY_TMAX);

        if (hit_rec.hit_anything)
        {
            // Material properties
            var material_type = i32(hit_rec.object_material.x);
            var fuzz = hit_rec.object_material.y;
            var ref_idx = hit_rec.object_material.z;
            var smoothness = hit_rec.object_material.w;

            var emitted = vec3f(0.0);

            // Generate random vectors for material functions
            var random_sphere = rng_next_vec3_in_unit_sphere(rng_state);

            var behaviour = material_behaviour(true, vec3f(0.0));

            // Dielectric 
            if (material_type < 0)
            {
                hit_rec.frontface = dot(ray_current.direction, hit_rec.normal) < 0.0;
                behaviour = dielectric(
                    hit_rec.normal,
                    ray_current.direction,
                    ref_idx,
                    hit_rec.frontface,
                    random_sphere,
                    fuzz,
                    rng_state
                );
            }
            // Metal (priority 3)
            else if (material_type == 1)
            {
                behaviour = metal(hit_rec.normal, ray_current.direction, fuzz, random_sphere);
            }
            // Lambertian (priority 4)
            else if (material_type == 0)
            {
                behaviour = lambertian(hit_rec.normal, 0.0, random_sphere, rng_state);
                light += color * hit_rec.object_color.rgb * hit_rec.object_material.w;

            }
            else
            {
                // Default to Lambertian if material type is unrecognized
                behaviour = lambertian(hit_rec.normal, 0.0, random_sphere, rng_state);
            }

            if (behaviour.scatter)
            {
                // Update accumulated color
                color = color * hit_rec.object_color.rgb;

                // Offset the hit point slightly along the normal to prevent self-intersection
                var offset_p = hit_rec.p + hit_rec.normal * RAY_TMIN;

                // Update the current ray
                ray_current = ray(offset_p, behaviour.direction);
            }
            else
            {
                // If the material doesn't scatter, add any emitted light and terminate
                break;
            }
        }
        else
        {
            // Ray didn't hit anything; apply environment color
            var env_color = envoriment_color(ray_current.direction, background_color1, background_color2);
            light = light + color * env_color;
            break;
        }
    }

    return light;
}

@compute @workgroup_size(THREAD_COUNT, THREAD_COUNT, 1)
fn render(@builtin(global_invocation_id) id : vec3u)
{
    var rez = uniforms[1];
    var time = u32(uniforms[0]);

    // Initialize RNG with pixel position, resolution, and frame
    var rng_state = init_rng(vec2u(id.x, id.y), vec2u(u32(rez)), time);

    // Get fragment coordinates
    var fragCoord = vec2f(f32(id.x), f32(id.y));
    var samples_per_pixel = i32(uniforms[4]);
    var color = vec3f(0.0);

    // Camera setup
    var lookfrom = vec3f(uniforms[7], uniforms[8], uniforms[9]);
    var lookat = vec3f(uniforms[23], uniforms[24], uniforms[25]);
    var cam = get_camera(lookfrom, lookat, vec3f(0.0, 1.0, 0.0), uniforms[10], 1.0, uniforms[6], uniforms[5]);

    // Loop for each sample per pixel
    for (var s = 0; s < samples_per_pixel; s = s + 1)
    {
        // Get UV coordinates with jitter for anti-aliasing
        var sample_uv = (fragCoord + sample_square(&rng_state)) / vec2f(rez);

        // Generate ray from camera through the pixel
        var r = get_ray(cam, sample_uv, &rng_state);

        // Trace the ray and accumulate color
        color = color + trace(r, &rng_state);  // Removed 'max_depth'
    }

    // Average the color over all samples
    color = color / f32(samples_per_pixel);

    // Apply gamma correction
    color = linear_to_gamma(color);
    var color_out = vec4f(color, 1.0);
    var map_fb = mapfb(id.xy, rez);

    // Accumulate the color if enabled
    var should_accumulate = uniforms[3];

    var new_color = should_accumulate * rtfb[map_fb] + color_out;
    rtfb[map_fb] = new_color;
    fb[map_fb] = new_color /rtfb[map_fb].w;

}