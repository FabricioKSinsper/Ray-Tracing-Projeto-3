fn hit_sphere(center: vec3f, radius: f32, r: ray, record: ptr<function, hit_record>, max: f32)
{
    var oc = r.origin - center;
    var a = dot(r.direction, r.direction);
    var half_b = dot(oc, r.direction);
    var c = dot(oc, oc) - radius * radius;
    var discriminant = half_b * half_b - a * c;

    if (discriminant < 0.0)
    {
        (*record).hit_anything = false;
        return;
    }

    var sqrt_discriminant = sqrt(discriminant);

    // Find the nearest root that lies within the acceptable range.
    var root = (-half_b - sqrt_discriminant) / a;
    if (root < RAY_TMIN || root > max)
    {
        root = (-half_b + sqrt_discriminant) / a;
        if (root < RAY_TMIN || root > max)
        {
            (*record).hit_anything = false;
            return;
        }
    }

    (*record).t = root;
    (*record).p = ray_at(r, root);
    var outward_normal = normalize((*record).p - center);

    // Determine if the ray is hitting the front face.
    (*record).frontface = dot(r.direction, outward_normal) < 0.0;
    if ((*record).frontface)
    {
        (*record).normal = outward_normal;
    }
    else
    {
        (*record).normal = -outward_normal;
    }
    (*record).hit_anything = true;
}

fn hit_quad(r: ray, Q: vec4f, u: vec4f, v: vec4f, record: ptr<function, hit_record>, max: f32)
{
  var n = cross(u.xyz, v.xyz);
  var normal = normalize(n);
  var D = dot(normal, Q.xyz);
  var w = n / dot(n.xyz, n.xyz);

  var denom = dot(normal, r.direction);
  if (abs(denom) < 0.0001)
  {
    record.hit_anything = false;
    return;
  }

  var t = (D - dot(normal, r.origin)) / denom;
  if (t < RAY_TMIN || t > max)
  {
    record.hit_anything = false;
    return;
  }

  var intersection = ray_at(r, t);
  var planar_hitpt_vector = intersection - Q.xyz;
  var alpha = dot(w, cross(planar_hitpt_vector, v.xyz));
  var beta = dot(w, cross(u.xyz, planar_hitpt_vector));

  if (alpha < 0.0 || alpha > 1.0 || beta < 0.0 || beta > 1.0)
  {
    record.hit_anything = false;
    return;
  }

  if (dot(normal, r.direction) > 0.0)
  {
    record.hit_anything = false;
    return;
  }

  record.t = t;
  record.p = intersection;
  record.normal = normal;
  record.hit_anything = true;
}

fn hit_triangle(r: ray, v0: vec3f, v1: vec3f, v2: vec3f, record: ptr<function, hit_record>, max: f32)
{
  var v1v0 = v1 - v0;
  var v2v0 = v2 - v0;
  var rov0 = r.origin - v0;

  var n = cross(v1v0, v2v0);
  var q = cross(rov0, r.direction);

  var d = 1.0 / dot(r.direction, n);

  var u = d * dot(-q, v2v0);
  var v = d * dot(q, v1v0);
  var t = d * dot(-n, rov0);

  if (u < 0.0 || u > 1.0 || v < 0.0 || (u + v) > 1.0)
  {
    record.hit_anything = false;
    return;
  }

  if (t < RAY_TMIN || t > max)
  {
    record.hit_anything = false;
    return;
  }

  record.t = t;
  record.p = ray_at(r, t);
  record.normal = normalize(n);
  record.hit_anything = true;
}

fn hit_box(r: ray, center: vec3f, rad: vec3f, record: ptr<function, hit_record>, t_max: f32)
{
  var m = 1.0 / r.direction;
  var n = m * (r.origin - center);
  var k = abs(m) * rad;

  var t1 = -n - k;
  var t2 = -n + k;

  var tN = max(max(t1.x, t1.y), t1.z);
  var tF = min(min(t2.x, t2.y), t2.z);

  if (tN > tF || tF < 0.0)
  {
    record.hit_anything = false;
    return;
  }

  var t = tN;
  if (t < RAY_TMIN || t > t_max)
  {
    record.hit_anything = false;
    return;
  }

  record.t = t;
  record.p = ray_at(r, t);
  record.normal = -sign(r.direction) * step(t1.yzx, t1.xyz) * step(t1.zxy, t1.xyz);
  record.hit_anything = true;

  return;
}

fn hit_cone(
    r: ray,
    apex: vec3f,
    base: vec3f,
    radius: f32,
    record: ptr<function, hit_record>,
    t_max: f32
) {
    // Vetor do ápice à base
    var axis = base - apex;
    var height = length(axis);
    var normalized_axis = normalize(axis);

    // Parâmetros do cone
    var cos_theta = radius / height;
    var cos_theta2 = cos_theta * cos_theta;

    // Origem do raio no espaço do cone
    var oc = r.origin - apex;

    // Coeficientes da equação quadrática
    var d_dot_a = dot(r.direction, normalized_axis);
    var oc_dot_a = dot(oc, normalized_axis);

    var a = dot(r.direction, r.direction) - cos_theta2 * d_dot_a * d_dot_a;
    var b = 2.0 * (dot(r.direction, oc) - cos_theta2 * d_dot_a * oc_dot_a);
    var c = dot(oc, oc) - cos_theta2 * oc_dot_a * oc_dot_a;

    var discriminant = b * b - 4.0 * a * c;

    // Verifica se há interseção
    if (discriminant < 0.0) {
        (*record).hit_anything = false;
        return;
    }

    var sqrt_discriminant = sqrt(discriminant);

    // Menor raiz
    var t = (-b - sqrt_discriminant) / (2.0 * a);
    if (t < RAY_TMIN || t > t_max) {
        t = (-b + sqrt_discriminant) / (2.0 * a);
        if (t < RAY_TMIN || t > t_max) {
            (*record).hit_anything = false;
            return;
        }
    }

    // Ponto de interseção
    var hit_point = ray_at(r, t);
    var height_at_hit = dot(hit_point - apex, normalized_axis);

    // Verifica se o ponto está dentro dos limites do cone
    if (height_at_hit < 0.0 || height_at_hit > height) {
        (*record).hit_anything = false;
        return;
    }

    // Calcula a normal
    var outward_normal = normalize(hit_point - (apex + height_at_hit * normalized_axis));
    (*record).normal = outward_normal;

    // Define o restante dos atributos do record
    (*record).t = t;
    (*record).p = hit_point;
    (*record).frontface = dot(r.direction, outward_normal) < 0.0;
    if (!(*record).frontface) {
        (*record).normal = -outward_normal;
    }
    (*record).hit_anything = true;
}