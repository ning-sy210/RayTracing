// FRAGMENT SHADER FOR SHADERTOY
// Run this at https://www.shadertoy.com/new
// See documentation at https://www.shadertoy.com/howto

// Your browser must support WebGL 2.0.
// Check your browser at https://webglreport.com/?v=2


//============================================================================
// Constants.
//============================================================================

const float PI = 3.1415926536;

const vec3 BACKGROUND_COLOR = vec3( 0.1, 0.2, 0.6 );

// Vertical field-of-view angle of camera. In degrees.
const float FOVY = 50.0;

// Use this for avoiding the "epsilon problem" or the shadow acne problem.
const float DEFAULT_TMIN = 10.0e-4;

// Use this for tmax for non-shadow ray intersection test.
const float DEFAULT_TMAX = 10.0e6;

// Equivalent to number of recursion levels (0 means ray-casting only).
// We are using iterations to replace recursions.
const int NUM_ITERATIONS = 10;

// Constants for the scene objects.
const int NUM_LIGHTS = 2;
const int NUM_MATERIALS = 5;
const int NUM_PLANES = 4;
const int NUM_SPHERES = 12;


//============================================================================
// Define new struct types.
//============================================================================
struct Ray_t {
    vec3 o;  // Ray Origin.
    vec3 d;  // Ray Direction. A unit vector.
};

struct Plane_t {
    // The plane equation is Ax + By + Cz + D = 0.
    float A, B, C, D;
    int materialID;
};

struct Sphere_t {
    vec3 center;
    float radius;
    int materialID;
};

struct Light_t {
    vec3 position;  // Point light 3D position.
    vec3 I_a;       // For Ambient.
    vec3 I_source;  // For Diffuse and Specular.
};

struct Material_t {
    vec3 k_a;   // Ambient coefficient.
    vec3 k_d;   // Diffuse coefficient.
    vec3 k_r;   // Reflected specular coefficient.
    vec3 k_rg;  // Global reflection coefficient.
    float n;    // The specular reflection exponent. Ranges from 0.0 to 128.0.
};

//----------------------------------------------------------------------------
// The lighting model used here is similar to that on Slides 8 and 12 of
// Lecture Topic 9 (Ray Tracing). Here it is computed as
//
//     I_local = SUM_OVER_ALL_LIGHTS {
//                   I_a * k_a +
//                   k_shadow * I_source * [ k_d * (N.L) + k_r * (R.V)^n ]
//               }
// and
//     I = I_local  +  k_rg * I_reflected
//----------------------------------------------------------------------------


//============================================================================
// Global scene data.
//============================================================================
Plane_t Plane[NUM_PLANES];
Sphere_t Sphere[NUM_SPHERES];
Light_t Light[NUM_LIGHTS];
Material_t Material[NUM_MATERIALS];



/////////////////////////////////////////////////////////////////////////////
// Initializes the scene.
/////////////////////////////////////////////////////////////////////////////
void InitScene()
{
    // Horizontal plane.
    Plane[0].A = 0.0;
    Plane[0].B = 1.0;
    Plane[0].C = 0.0;
    Plane[0].D = 0.0;
    Plane[0].materialID = 0;

    // Vertical plane (x-y).
    Plane[1].A = 0.0;
    Plane[1].B = 0.0;
    Plane[1].C = 1.0;
    Plane[1].D = 3.5;
    Plane[1].materialID = 0;
    
    // Vertical plane (y-z).
    Plane[2].A = 1.0;
    Plane[2].B = 0.0;
    Plane[2].C = 0.0;
    Plane[2].D = 3.5;
    Plane[2].materialID = 0;
    
    // Vertical plane (y-z).
    Plane[3].A = 1.0;
    Plane[3].B = 0.0;
    Plane[3].C = 0.0;
    Plane[3].D = -3.5;
    Plane[3].materialID = 0;

    // Rotating sphere 1
    Sphere[0].center = vec3( sin(4.0 * iTime), cos(4.0 * iTime) + 1.4, cos(4.0 * iTime) );
    Sphere[0].radius = 0.13;
    Sphere[0].materialID = 1;
    
    // Rotating sphere 2
    Sphere[1].center = vec3( -sin(4.0 * iTime + PI / 6.0), -sin(4.0 * iTime + PI / 6.0) + 1.4, -cos(4.0 * iTime + PI / 6.0) );
    Sphere[1].radius = 0.13;
    Sphere[1].materialID = 1;
    
    // Rotating sphere 3
    Sphere[2].center = vec3( sin(4.0 * iTime + PI / 3.0), 1.4, cos(4.0 * iTime + PI / 3.0) );
    Sphere[2].radius = 0.13;
    Sphere[2].materialID = 1;
    
    // Stationary Sphere 1
    Sphere[3].center = vec3( 0, 1.4, 0 );
    Sphere[3].radius = 0.12;
    Sphere[3].materialID = 3;
    
    // Stationary Sphere 2
    Sphere[4].center = vec3( 0.1, 1.2, 0.1 );
    Sphere[4].radius = 0.12;
    Sphere[4].materialID = 2;
    
    // Stationary Sphere 3
    Sphere[5].center = vec3( -0.1, 1.4, -0.1 );
    Sphere[5].radius = 0.12;
    Sphere[5].materialID = 2;

    // Stationary Sphere 4
    Sphere[6].center = vec3( 0.14, 1.4, -0.15 );
    Sphere[6].radius = 0.12;
    Sphere[6].materialID = 3;

    // Stationary Sphere 5
    Sphere[7].center = vec3( 0.18, 1.4, 0.18);
    Sphere[7].radius = 0.12;
    Sphere[7].materialID = 2;

    // Stationary Sphere 6
    Sphere[8].center = vec3( -0.18, 1.4, 0.18);
    Sphere[8].radius = 0.12;
    Sphere[8].materialID = 2;

    // Stationary Sphere 6
    Sphere[9].center = vec3( 0.04, 1.24, 0.34);
    Sphere[9].radius = 0.12;
    Sphere[9].materialID = 3;

    // Stationary Sphere 8
    Sphere[10].center = vec3( 0.04, 1.5, 0.34);
    Sphere[10].radius = 0.12;
    Sphere[10].materialID = 3;
    
    // Stationary Sphere 9
    Sphere[11].center = vec3( 0.05, 1.62, 0.0);
    Sphere[11].radius = 0.12;
    Sphere[11].materialID = 2;

    // Silver material.
    Material[0].k_d = vec3( 0.5, 0.5, 0.5 );
    Material[0].k_a = 0.2 * Material[0].k_d;
    Material[0].k_r = 2.0 * Material[0].k_d;
    Material[0].k_rg = 0.5 * Material[0].k_r;
    Material[0].n = 64.0;

    // Black material.
    Material[1].k_d = vec3( 0.34, 0.34, 0.34 );
    Material[1].k_a = 0.2 * Material[1].k_d;
    Material[1].k_r = 0.5 * Material[1].k_d;
    Material[1].k_rg = 0.5 * Material[1].k_r;
    Material[1].n = 64.0;

    // Red material.
    Material[2].k_d = vec3( 0.83, 0.29, 0.27 );
    Material[2].k_a = 0.2 * Material[2].k_d;
    Material[2].k_r = vec3( 0.4, 0.4, 0.4 );
    Material[2].k_rg = 0.5 * Material[2].k_r;
    Material[2].n = 64.0;
    
    // Blue material.
    Material[3].k_d = vec3( 0.1, 0.82, 0.96 );
    Material[3].k_a = 0.2 * Material[2].k_d;
    Material[3].k_r = vec3( 0.4, 0.4, 0.4 );
    Material[3].k_rg = 0.5 * Material[2].k_r;
    Material[3].n = 128.0;

    // Light 0.
    Light[0].position = vec3( 3.0 * sin(5.0 * iTime), 12, 3.0 * cos(5.0 * iTime) );
    Light[0].I_a = vec3( 0.5, 0.5, 0.5 );
    Light[0].I_source = vec3( 0.5, cos(iTime), sin(iTime) );
    
    // Light 1.
    Light[1].position = vec3( 3.0 * cos(3.0 * iTime), 8.0, 3.0 * sin(3.0 * iTime) );
    Light[1].I_a = vec3( sin(iTime), cos(iTime), 0.8 );
    Light[1].I_source = vec3( 0.4, 0.4, 0.4 );
}



/////////////////////////////////////////////////////////////////////////////
// Returns a random number between 0 and 1.
//
// This pseudorandom number generator is based on the 32-bit combined LFSR
// generator proposed in the paper "Tables of Maximally-Equidistributed
// Combined LFSR Generators" by Pierre L'Ecuyer.
// (http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.43.3639)
/////////////////////////////////////////////////////////////////////////////

// VERY IMPORTANT: The initial seeds rand_z1, rand_z2, rand_z3, rand_z4
// must be larger than 1, 7, 15, and 127 respectively.
const uint CONST_RAND_SEED = 987654321U;
uint rand_z1 = uint(CONST_RAND_SEED + 2U);
uint rand_z2 = uint(CONST_RAND_SEED + 8U);
uint rand_z3 = uint(CONST_RAND_SEED + 16U);
uint rand_z4 = uint(CONST_RAND_SEED + 128U);

float rand(void)
{
    uint b  = ((rand_z1 << 6) ^ rand_z1) >> 13;
    rand_z1 = ((rand_z1 & 4294967294U) << 18) ^ b;
    b       = ((rand_z2 << 2) ^ rand_z2) >> 27;
    rand_z2 = ((rand_z2 & 4294967288U) << 2) ^ b;
    b       = ((rand_z3 << 13) ^ rand_z3) >> 21;
    rand_z3 = ((rand_z3 & 4294967280U) << 7) ^ b;
    b       = ((rand_z4 << 3) ^ rand_z4) >> 12;
    rand_z4 = ((rand_z4 & 4294967168U) << 13) ^ b;
    return float(rand_z1 ^ rand_z2 ^ rand_z3 ^ rand_z4) * 2.3283064365386963e-10;
}



/////////////////////////////////////////////////////////////////////////////
// Computes intersection between a plane and a ray.
// Returns true if there is an intersection where the ray parameter t is
// between tmin and tmax, otherwise returns false.
// If there is such an intersection, outputs the value of t, the position
// of the intersection (hitPos) and the normal vector at the intersection
// (hitNormal).
/////////////////////////////////////////////////////////////////////////////
bool IntersectPlane( in Plane_t pln, in Ray_t ray,
                     in float tmin, in float tmax,
                     out float t, out vec3 hitPos, out vec3 hitNormal )
{
    vec3 N = vec3( pln.A, pln.B, pln.C );
    float NRd = dot( N, ray.d );
    float NRo = dot( N, ray.o );
    float t0 = (-pln.D - NRo) / NRd;
    if ( t0 < tmin || t0 > tmax ) return false;

    // We have a hit -- output results.
    t = t0;
    hitPos = ray.o + t0 * ray.d;
    hitNormal = normalize( N );
    return true;
}



/////////////////////////////////////////////////////////////////////////////
// Computes intersection between a plane and a ray.
// Returns true if there is an intersection where the ray parameter t is
// between tmin and tmax, otherwise returns false.
/////////////////////////////////////////////////////////////////////////////
bool IntersectPlane( in Plane_t pln, in Ray_t ray,
                     in float tmin, in float tmax )
{
    vec3 N = vec3( pln.A, pln.B, pln.C );
    float NRd = dot( N, ray.d );
    float NRo = dot( N, ray.o );
    float t0 = (-pln.D - NRo) / NRd;
    if ( t0 < tmin || t0 > tmax ) return false;
    return true;
}



/////////////////////////////////////////////////////////////////////////////
// Computes intersection between a sphere and a ray.
// Returns true if there is an intersection where the ray parameter t is
// between tmin and tmax, otherwise returns false.
// If there is one or two such intersections, outputs the value of the
// smaller t, the position of the intersection (hitPos) and the normal
// vector at the intersection (hitNormal).
/////////////////////////////////////////////////////////////////////////////
bool IntersectSphere( in Sphere_t sph, in Ray_t ray,
                      in float tmin, in float tmax,
                      out float t, out vec3 hitPos, out vec3 hitNormal )
{
    ///////////////////////////////////
    // TASK 1: WRITE YOUR CODE HERE. //
    ///////////////////////////////////

    // need to account for the fact that the sphere is not at origin
    float a = dot(ray.d, ray.d);
    float b = 2.0 * (dot(ray.o, ray.d) - dot(ray.d, sph.center));
    float c = dot(ray.o, ray.o) - sph.radius * sph.radius + dot(sph.center, sph.center) - 2.0 * dot(ray.o, sph.center);

    float d = b * b - 4.0 * a * c;
    float t0;

    if (d > 0.0) {                                  // 2 intersection points
        float t1 = (-b + sqrt(d)) / (2.0 * a);
        float t2 = (-b - sqrt(d)) / (2.0 * a);

        if (max(t1, t2) < tmin || min(t1, t2) > tmax) {
            return false;
        } else {
            t0 = min(t1, t2);
        }
    } else if (d == 0.0) {                          // only 1 intersection point 
        if (t0 < tmin || t0 > tmax) {
            return false;
        }

        t0 = (-b + sqrt(d)) / (2.0 * a);
    } else {                                        // no real intersection points
        return false;
    }

    t = t0;
    hitPos = ray.o + t0 * ray.d;
    hitNormal = normalize(hitPos - sph.center);

    return true;
}



/////////////////////////////////////////////////////////////////////////////
// Computes intersection between a sphere and a ray.
// Returns true if there is an intersection where the ray parameter t is
// between tmin and tmax, otherwise returns false.
/////////////////////////////////////////////////////////////////////////////
bool IntersectSphere( in Sphere_t sph, in Ray_t ray,
                      in float tmin, in float tmax )
{
    ///////////////////////////////////
    // TASK 1: WRITE YOUR CODE HERE. //
    ///////////////////////////////////

    float a = dot(ray.d, ray.d);
    float b = 2.0 * (dot(ray.o, ray.d) - dot(ray.d, sph.center));
    float c = dot(ray.o, ray.o) - sph.radius * sph.radius + dot(sph.center, sph.center) - 2.0 * dot(ray.o, sph.center);

    float d = b * b - 4.0 * a * c;

    if (d > 0.0) {                                                      // 2 intersection points
        float t1 = (-b + sqrt(d)) / (2.0 * a);
        float t2 = (-b - sqrt(d)) / (2.0 * a);

        return !((t1 < tmin || t1 > tmax) && (t2 < tmin || t2 > tmax));
    } else if (d == 0.0) {                                              // only 1 intersection point 
        float temp = (-b + sqrt(d)) / (2.0 * a);

        return !(temp < tmin || temp > tmax);
    } else {                                                            // no real intersection points
        return false;
    }
}



/////////////////////////////////////////////////////////////////////////////
// Computes (I_a * k_a) + k_shadow * I_source * [ k_d * (N.L) + k_r * (R.V)^n ].
// Input vectors L, N and V are pointing AWAY from surface point.
// Assume all vectors L, N and V are unit vectors.
/////////////////////////////////////////////////////////////////////////////
vec3 PhongLighting( in vec3 L, in vec3 N, in vec3 V, in bool inShadow,
                    in Material_t mat, in Light_t light )
{
    if ( inShadow ) {
        return light.I_a * mat.k_a;
    }
    else {
        vec3 R = reflect( -L, N );
        float N_dot_L = max( 0.0, dot( N, L ) );
        float R_dot_V = max( 0.0, dot( R, V ) );
        float R_dot_V_pow_n = ( R_dot_V == 0.0 )? 0.0 : pow( R_dot_V, mat.n );

        return light.I_a * mat.k_a +
               light.I_source * (mat.k_d * N_dot_L + mat.k_r * R_dot_V_pow_n);
    }
}



/////////////////////////////////////////////////////////////////////////////
// Casts a ray into the scene and returns color computed at the nearest
// intersection point. The color is the sum of light from all light sources,
// each computed using Phong Lighting Model, with consideration of
// whether the interesection point is being shadowed from the light.
// If there is no interesection, returns the background color, and outputs
// hasHit as false.
// If there is intersection, returns the computed color, and outputs
// hasHit as true, the 3D position of the intersection (hitPos), the
// normal vector at the intersection (hitNormal), and the k_rg value
// of the material of the intersected object.
/////////////////////////////////////////////////////////////////////////////
vec3 CastRay( in Ray_t ray,
              out bool hasHit, out vec3 hitPos,
              out vec3 hitNormal, out vec3 k_rg )
{
    // Find whether and where the ray hits some object.
    // Take the nearest hit point.

    bool hasHitSomething = false;
    float nearest_t = DEFAULT_TMAX;   // The ray parameter t at the nearest hit point.
    vec3 nearest_hitPos;              // 3D position of the nearest hit point.
    vec3 nearest_hitNormal;           // Normal vector at the nearest hit point.
    int nearest_hitMatID;             // MaterialID of the object at the nearest hit point.

    float temp_t;
    vec3 temp_hitPos;
    vec3 temp_hitNormal;
    bool temp_hasHit;

    ///////////////////////////////////////////////////////////////////////////
    // * Intersect input ray with all the planes and spheres, and record the
    //   front-most (nearest) intersection.
    // * If there is intersection, need to record hasHitSomething,
    //   nearest_t, nearest_hitPos, nearest_hitNormal, nearest_hitMatID.
    ///////////////////////////////////////////////////////////////////////////

    // test for intersection with all planes
    for (int i = 0; i < NUM_PLANES; i++) {
        temp_hasHit = IntersectPlane(Plane[i], ray, DEFAULT_TMIN, DEFAULT_TMAX, temp_t, temp_hitPos, temp_hitNormal);

        if (temp_hasHit) {
            if (!hasHitSomething) {
                hasHitSomething = true;
            }

            if (temp_t < nearest_t) {
                nearest_t = temp_t;
                nearest_hitPos = temp_hitPos;
                nearest_hitNormal = temp_hitNormal;
                nearest_hitMatID = Plane[i].materialID;
            }
        }
    }

    // test for intersection with all spheres
    for (int i = 0; i < NUM_SPHERES; i++) {
        temp_hasHit = IntersectSphere(Sphere[i], ray, DEFAULT_TMIN, DEFAULT_TMAX, temp_t, temp_hitPos, temp_hitNormal);

        if (temp_hasHit) {
            if (!hasHitSomething) {
                hasHitSomething = true;
            }

            if (temp_t < nearest_t) {
                nearest_t = temp_t;
                nearest_hitPos = temp_hitPos;
                nearest_hitNormal = temp_hitNormal;
                nearest_hitMatID = Sphere[i].materialID;
            }
        }
    }

    // One of the output results.
    hasHit = hasHitSomething;
    if ( !hasHitSomething ) return BACKGROUND_COLOR;

    vec3 I_local = vec3( 0.0 );  // Result color will be accumulated here.

    ///////////////////////////////////////////////////////////////////////////
    // * Accumulate lighting from each light source on the nearest hit point.
    //   They are all accumulated into I_local.
    // * For each light source, make a shadow ray, and check if the shadow ray
    //   intersects any of the objects (the planes and spheres) between the
    //   nearest hit point and the light source.
    // * Then, call PhongLighting() to compute lighting for this light source.
    ///////////////////////////////////////////////////////////////////////////

    for (int i = 0; i < NUM_LIGHTS; i++) {
        vec3 cam_pos = vec3( 2.5, 1.0, 2.5 );
        
        // make shadow ray for the intersection point towards each light source
        Ray_t shadow_ray;
        shadow_ray.o = nearest_hitPos;
        shadow_ray.d = normalize(Light[i].position - nearest_hitPos);

        vec3 L = shadow_ray.d;
        vec3 N = nearest_hitNormal;
        vec3 V = normalize(cam_pos - nearest_hitPos);

        float shadow_ray_max_t = distance(Light[i].position, nearest_hitPos);

        // check intersection with planes
        for (int j = 0; j < NUM_PLANES; j++) {
            temp_hasHit = IntersectPlane(Plane[j], shadow_ray, DEFAULT_TMIN, shadow_ray_max_t);
            
            // no need to look for further intersections once we know that the intersection point is in shadow
            if (temp_hasHit) {
                break;
            }
        }

        // check intersection with spheres only if the shadow ray has not intersected with any planes yet
        if (!temp_hasHit) {
            for (int k = 0; k < NUM_SPHERES; k++) {
                temp_hasHit = IntersectSphere(Sphere[k], shadow_ray, DEFAULT_TMIN, shadow_ray_max_t);

                // no need to look for further intersections once we know that the intersection point is in shadow
                if (temp_hasHit) {
                    break;
                }
            }    
        }
        
        I_local += PhongLighting(L, N, V, temp_hasHit, Material[nearest_hitMatID], Light[i]);
    }

    // Populate output results.
    hitPos = nearest_hitPos;
    hitNormal = nearest_hitNormal;
    k_rg = Material[nearest_hitMatID].k_rg;

    return I_local;
}



/////////////////////////////////////////////////////////////////////////////
// Execution of fragment shader starts here.
// 1. Initializes the scene.
// 2. Compute a primary ray for the current pixel (fragment).
// 3. Trace ray into the scene with NUM_ITERATIONS recursion levels.
/////////////////////////////////////////////////////////////////////////////
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Initialize random number generator before the first call to rand().
    uint RAND_SEED = uint( (mod(iTime*100.0, 100.0) + 101.01) *
                           (fragCoord.x + 17.0) * (fragCoord.y + 23.0) );
    rand_z1 = uint(RAND_SEED + 2U);
    rand_z2 = uint(RAND_SEED + 8U);
    rand_z3 = uint(RAND_SEED + 16U);
    rand_z4 = uint(RAND_SEED + 128U);


    InitScene();


    // Camera position and orientation in world space.
    vec3 cam_pos = vec3( 3.0, 1.4, 5.0 );
    vec3 cam_lookat = vec3( 0.0, 1.4, -1.0 );
    vec3 cam_up_vec = vec3( 0.0, 1.0, 0.0 );

    // Camera coordinate frame in world space.
    vec3 cam_z_axis = normalize( cam_pos - cam_lookat );
    vec3 cam_x_axis = normalize( cross(cam_up_vec, cam_z_axis) );
    vec3 cam_y_axis = normalize( cross(cam_z_axis, cam_x_axis));

    // Vertical field-of-view angle of camera. In radians.
    float cam_FOVY = FOVY * PI / 180.0;

    // Perpendicular distance of the image rectangle from the camera.
    // If implementing depth-of-field, the plane of the image rectangle
    // is the plane of focus.
    float image_dist = distance(cam_pos, vec3(0.0, 0.7, 0.0));

    float image_height = 2.0 * image_dist * tan(cam_FOVY / 2.0);
    float image_width = image_height * iResolution.x / iResolution.y;
    float pixel_width = image_width / iResolution.x;

    // Image rectangle origin (bottom-leftmost corner) position in camera space.
    vec3 image_origin = vec3(-image_width/2.0, -image_height/2.0, -image_dist);


    /////////////////////////////////////////////////////////////////////////
    // * Trace multiple (SPP) random primary rays per pixel to produce
    //   depth-of-field effect and for image anti-aliasing (reduce jaggies).
    // * Each primary ray starts from a random position on the lens and
    //   points towards a random position inside the current pixel.
    // * The lens is assumed to have a square-shaped aperture of size
    //   aperture_width x aperture_width. The lens is centered at the
    //   the origin of the camera frame, and parallel to the x-y plane
    //   of the camera frame.
    // * The final color of the current pixel is the average color of
    //   all the primary rays.
    /////////////////////////////////////////////////////////////////////////

    //=======================================================================
    // These constants are used for distribution ray tracing to produce
    // depth-of-field effect and for image anti-aliasing (reduce jaggies).
    //=======================================================================
    // Number of samples (random primary rays) per pixel.
    const int SPP = 32;

    // Lens aperture width. Assume square aperture.
    const float aperture_width = 0.3;
    //=======================================================================

    // Current pixel 3D position in camera space.
    vec3 pixel_pos = image_origin +
                     vec3(pixel_width * fragCoord, 0);

    // Create primary ray.
    Ray_t pRay;
    pRay.o = cam_pos;
    pRay.d = normalize( cam_pos + pixel_pos.x * cam_x_axis  +
                                  pixel_pos.y * cam_y_axis  +
                                  pixel_pos.z * cam_z_axis - pRay.o);

    // Start Ray Tracing.
    // Use iterations to emulate the recursion.

    vec3 I_result = vec3( 0.0 );
    vec3 compounded_k_rg = vec3( 1.0 );
    Ray_t nextRay = pRay;

    for ( int level = 0; level <= NUM_ITERATIONS; level++ )
    {
        bool hasHit;
        vec3 hitPos, hitNormal, k_rg;

        vec3 I_local = CastRay( nextRay, hasHit, hitPos, hitNormal, k_rg );

        I_result += compounded_k_rg * I_local;

        if ( !hasHit ) break;

        compounded_k_rg *= k_rg;

        nextRay = Ray_t( hitPos, normalize( reflect(nextRay.d, hitNormal) ) );
    }

    fragColor = vec4( I_result, 1.0 );
}
