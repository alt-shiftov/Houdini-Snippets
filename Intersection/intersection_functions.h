/////// LINE - LINE
// https://paulbourke.net/geometry/pointlineplane/lineline.c

// intersection between lines - is just finding minimal distance between two lines
int line_line_intersect(vector p1, p2, p3, p4, pa, pb; float parameter[]){
    
    float EPS = 1e-6;
    vector p13,p43,p21;
    
    // Vectors (p13, p43, p21) representing differences between points
    p13 = set(p1.x - p3.x, p1.y - p3.y, p1.z - p3.z);
    p43 = set(p4.x - p3.x, p4.y - p3.y, p4.z - p3.z);
    
    // Check for degenerate cases (segments are points)
    if (abs(p43.x) < EPS && abs(p43.y) < EPS && abs(p43.z) < EPS)
      return 0;
    
    p21 = set(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
    
    if (abs(p21.x) < EPS && abs(p21.y) < EPS && abs(p21.z) < EPS)
      return 0;
      
    // Calculate dot products for intersection check
    float d1343 = p13.x * p43.x + p13.y * p43.y + p13.z * p43.z;
    float d4321 = p43.x * p21.x + p43.y * p21.y + p43.z * p21.z;
    float d1321 = p13.x * p21.x + p13.y * p21.y + p13.z * p21.z;
    float d4343 = p43.x * p43.x + p43.y * p43.y + p43.z * p43.z;
    float d2121 = p21.x * p21.x + p21.y * p21.y + p21.z * p21.z;
    
    
    // Calculate denominator for intersection parameter
    float denom = d2121 * d4343 - d4321 * d4321;
    if (abs(denom) < EPS)
      return 0;
      
     // Calculate numerator and intersection parameters
    float numer = d1343 * d4321 - d1321 * d4343;
    float mua = numer / denom;
    float mub = (d1343 + d4321 * (mua)) / d4343;
    
    // Store parametric values in the parameter array
    resize(parameter, 2);
    parameter[0] = mua;
    parameter[1] = mub;
    
    
    // Calculate intersection points
    pa = set(p1.x + mua * p21.x, p1.y + mua * p21.y, p1.z + mua * p21.z);
    pb = set(p3.x + mub * p43.x, p3.y + mub * p43.y, p3.z + mub * p43.z);
    
    
    
    return 1;
}

////////////////////////////////////////////
/////// SEGMENT - SEGMENT
// https://paulbourke.net/geometry/pointlineplane/lineline.c

// The intersect_segment_segment() function is a variation of line_line_intersect(), 
// but it additionally restricts the intersection point to lie 
// within the parametric bounds of both segments (0 <= mua, mub <= 1). 
// By adjusting these bounds, we can adapt the function to handle intersections 
// between rays, lines, and segments in various combinations

int intersect_segment_segment(vector p1, p2, p3, p4, pa, pb; float parameter[]){
    vector p13, p43, p21;
    float EPS = 1e-6;
    
    // Vectors (p13, p43, p21) representing differences between points
    p13 = set(p1.x - p3.x, p1.y - p3.y, p1.z - p3.z);
    p43 = set(p4.x - p3.x, p4.y - p3.y, p4.z - p3.z);
    
    // Check for degenerate cases (segments are points)
    if (abs(p43.x) < EPS && abs(p43.y) < EPS && abs(p43.z) < EPS)
      return 0;
    
    p21 = set(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
    
    if (abs(p21.x) < EPS && abs(p21.y) < EPS && abs(p21.z) < EPS)
      return 0;
      
    // Calculate dot products for intersection check
    float d1343 = p13.x * p43.x + p13.y * p43.y + p13.z * p43.z;
    float d4321 = p43.x * p21.x + p43.y * p21.y + p43.z * p21.z;
    float d1321 = p13.x * p21.x + p13.y * p21.y + p13.z * p21.z;
    float d4343 = p43.x * p43.x + p43.y * p43.y + p43.z * p43.z;
    float d2121 = p21.x * p21.x + p21.y * p21.y + p21.z * p21.z;
    
    // Calculate denominator for intersection parameter
    float denom = d2121 * d4343 - d4321 * d4321;
    if (abs(denom) < EPS)
      return 0;
      
    // Calculate numerator and intersection parameters
    float numer = d1343 * d4321 - d1321 * d4343;
    float mua = numer / denom;
    float mub = (d1343 + d4321 * (mua)) / d4343;
    
    // Check if intersection points are within segment bounds
    if (mua < 0.0 || mua > 1.0 || mub < 0.0 || mub > 1.0) {
      pa = set(0,0,0); // Set intersection points to zero if outside bounds
      pb = set(0,0,0);
      return 0;
    }
    
    // Calculate intersection points
    pa = set(p1.x + mua * p21.x, p1.y + mua * p21.y, p1.z + mua * p21.z);
    pb = set(p3.x + mub * p43.x, p3.y + mub * p43.y, p3.z + mub * p43.z);
    
    // Store parametric values in the parameter array
    resize(parameter, 2);
    parameter[0] = mua;
    parameter[1] = mub;
    
    // Check if intersection points are close enough
    if (distance(pa, pb) < 0.001)
        return 1; // Intersection found
    else
        return 0; // No intersection
}


////////////////////////////////////////
///// SEGMENT - CURVE

// Curve - is multiple segments.

int boxes_intersect2d(vector A_min, A_max, B_min, B_max)
{
    // Check overlap on x-axis
    if (A_max.x < B_min.x || B_max.x < A_min.x) {
        return 0; // No overlap
    }

    // Check overlap on y-axis
    if (A_max.y < B_min.y || B_max.y < A_min.y) {
        return 0; // No overlap
    }

    // // Check overlap on z-axis
    // if (A_max.z < B_min.z || B_max.z < A_min.z) {
    //     return 0; // No overlap
    // }

    return 1; // Intersection
}

// сделать curve через массив позиций
int intersect_segment_curve(int geo; int prim; vector p0, p1; vector positions[], uvw[]; float parameter[]){
    // Finds intersections between a segment defined by points p0 and p1, and a curve defined by geometry geo and primitive prim.
    // Returns the number of intersections found.
    // Populates positions[] with intersection points, uvw[] with curve parameterization at intersections, and parameter[] with segment parameters.

    int pts[] = primpoints(geo, prim); // get curve points
    float eps = 1e-6;

    // vector A_max = max(p0, p1);
    // vector A_min = min(p0, p1);

    // Optimize bounding box calculation using quaternion for alignment
    vector4 quat = dihedral(p1-p0, set(1,0,0)); // Quaternion for aligning segment by Z-axis
    
    vector A_min = qrotate(quat, p0);
    vector A_max = qrotate(quat, p1);

    // Iterate over curve segments
    for(int i = 1; i < len(pts); i++){
        vector q0 = point(geo, "P", pts[i-1]);
        vector q1 = point(geo, "P", pts[i]);
       
        // Orient curve segment bounds using the same quaternion
        vector q0_rot = qrotate(quat, q0);
        vector q1_rot = qrotate(quat, q1);

        // Calculate bounding box for the current curve segment
        vector B_max = max(q0_rot, q1_rot);
        vector B_min = min(q0_rot, q1_rot);

        // Early rejection: check if bounding boxes intersect
        int isinside = boxes_intersect2d(A_min, A_max, B_min, B_max);
        if (!isinside) continue; // If bounding boxes don't intersect, skip to next segment
        
        vector P; float params[];
        // Check for intersection between the input segment and the current curve segment
        int hasinters = intersect_segment_segment(p0, p1, q0, q1, P, params);
        if (hasinters){
            // Append intersection point to output array
            append(positions, P);
            
            // Calculate U coordinate based on curve parameterization
            float u0 = (i-1) / (len(pts)*1.0-1);
            float u1 = i / (len(pts)*1.0-1);

            // Using parameter for lerp between u0 and u1
            float u = lerp(u0, u1, params[1]);

            append(uvw, set(u, 0, 0)); // uvw - curve intersection barycentric coords
            append(parameter, params[0]); // parameter - segment intersection
        }
        
    }
    return len(positions);
    
}



////////////////////////////////////////////
////// LINE - PLANE (PLANE is described by normal and point)

function vector line_plane_intersection(vector planeNormal, planePoint, rayDirection, rayPoint; vector interspos){
    // Calculates the intersection point of a line (defined by rayDirection and rayPoint)
    // with a plane (defined by planeNormal and planePoint).

    float epsilon = 1e-6;
    
    // Calculate the dot product of the plane normal and the ray direction
    float ndotu = dot(planeNormal, rayDirection);
    
    // Check for parallel or coincident line and plane
    if (abs(ndotu) < epsilon){
            // no intersection or line is within plane
            return 0;
    }
    
    // Calculate the vector from the plane point to the ray point
    vector w = rayPoint - planePoint;
    
    // Calculate the distance from the ray point to the intersection point
    float si = dot(-planeNormal, w) / ndotu;
    
    // Calculate the intersection point
    interspos = w + si * rayDirection + planePoint;
    return 1; // Intersection point found
}


////// LINE - PLANE (PLANE is described by three points)
// https://paulbourke.net/geometry/pointlineplane/pointofintersection.txt

vector CrossProduct(vector p0, p1, p2){ 
    // Calculates the normal vector of a plane defined by three points.
    // Calculate two vectors in the plane
    vector A = set(p1.x - p0.x, p1.y - p0.y, p1.z - p0.z);
    vector B = set(p2.x - p0.x, p2.y - p0.y, p2.z - p0.z);
    
     // Calculate the cross product of the two vectors to get the normal vector
    return cross(A, B);
}


int line_plane_intersection(vector La, Lb, P0, P1, P2, intersection)
{
    // thanks to Paul Bourke and Bryan Hanson
    // Finds the intersection point of a line defined by points La and Lb
    // with a plane defined by points P0, P1, and P2.

    // Calculate the normal vector of the plane
    vector N = CrossProduct(P0, P1, P2); 
    
    // Calculate the denominator for the line parameter
    float n = dot(N, set(P2.x - La.x, P2.y - La.y, P2.z - La.z));
    vector LbMinusLa = set(Lb.x - La.x, Lb.y - La.y, Lb.z - La.z);
    float d = dot(N, LbMinusLa);
    
    // Check for parallel or coincident line and plane
    if (d == 0)
        return -1; //Line is parallel to or in the plane
    
    // Calculate the line parameter
    float u = n / d;
    
    // Check if the intersection point is within the line segment (optional)
    // if ((u >= 0.0) && (u <= 1.0)) {
    //     // Plane is between the two points
    // }
    
    // Calculate the intersection point
    intersection = set(La.x + u * LbMinusLa.x, La.y + u * LbMinusLa.y, La.z + u * LbMinusLa.z);
    
    return 1; // Intersection point found
}


////////////////////////////////////////////
////// PLANE - PLANE Intersection. Planes is described by point and normal
// https://paulbourke.net/geometry/pointlineplane/

function int planes_intersection(vector n1, p1, n2, p2, point1, point2)
{
    float eps = 1e-6;

    // Check for parallel planes
    if (length(cross(n1, n2)) < eps){
        // planes are parallel
        return 0;
    }
    
    // Calculate determinant for plane equation coefficients
    float det = dot(n1, n1)*dot(n2, n2) - pow((dot(n1, n2)), 2);
    
    // Calculate plane distances from origin
    float d1 = dot(n1, p1);
    float d2 = dot(n2, p2);
    
    // Calculate coefficients for intersection line
    float c1 = (dot(d1*n2, n2) - dot(d2*n1, n2)) / det;
    float c2 = (dot(d2*n1, n1) - dot(d1*n1, n2)) / det;
    
    // Calculate a point on the intersection line
    float u = 0;
    vector ptpos = c1*n1 + c2*n2 + u*n1*n2;
    
    // Calculate direction vector of the intersection line
    vector dir = cross(n1, n2);
    dir = normalize(dir);
    
    // Calculate two points on the intersection line
    float distance = 10.0; // Adjust as needed
    point1 = ptpos + distance * normalize(dir);
    point2 = ptpos - distance * normalize(dir);
    
    return 1; // Planes intersect
    
}


////////////////////////////////////////////
///// SEGMENT - TRIANGLE Intersection

// A,B,C - triangle points
// P0, P1 - segment points
// P - result position

int intersect_segment_triangle(vector A, B, C, P0, P1, P; float parameter; vector uvw){
    
    float eps = 1e-6;
    
    // Edge vectors for triangle
    vector AB = B - A;
    vector AC = C - A;
    
    // Plane normal
    vector normal = cross(AB, AC);
    
    // Line direction
    vector dir = P1 - P0;
    
    // Check if line and plane are parallel
    float dot_nd = dot(normal, dir);
    if (abs(dot_nd) < eps) {
        return 0; // No intersection
    }
    
    // Compute distance from P0 to plane
    float d = dot(normal, A - P0) / dot_nd;
    
     // Check if intersection point is on the line segment
    if (d < 0.0 || d > 1.0) {
        return 0; // No intersection on the line segment
    }
    
    // Intersection point on line
    P = P0 + d * dir;
    
    // Barycentric coordinates
    vector v0 = B - A;
    vector v1 = C - A;
    vector v2 = P - A;
    
    float dot00 = dot(v0, v0);
    float dot01 = dot(v0, v1);
    float dot02 = dot(v0, v2);
    float dot11 = dot(v1, v1);
    float dot12 = dot(v1, v2);
    
    // Compute barycentric coordinates (u, v, w)
    float denom = dot00 * dot11 - dot01 * dot01;
    if (abs(denom) < eps) {
        return 0; // No intersection or degenerate triangle
    }
    
    float u = (dot11 * dot02 - dot01 * dot12) / denom;
    if (u < 0.0 || u > 1.0) {
        return 0; // Intersection outside triangle
    }
    
    float v = (dot00 * dot12 - dot01 * dot02) / denom;
    if (v < 0.0 || u + v > 1.0) {
        return 0; // Intersection outside triangle
    }
    
    parameter = d; // Parametric position (0 <= d <= 1) on the segment
    uvw = set(u, v, 0); // Barycentric coords of intersection on the triangle
    
    // Intersection point is within the triangle
    return 1;
}


///////////////////////////////////////
////// TRIANGLE - TRIANGLE INTERSECTION

function void triangles_intersection(vector A0, B0, C0, A1, B1, C1; vector positions[], uvw0[], uvw1[]){
    // Function to find the intersection points between two triangles.
    // It iterates through all edges of both triangles, checking for intersections
    // with the opposite triangle using intersect_segment_triangle().

    vector pos, uvw;
    float parameter;
    
    int i = 0;

    // how vertices is described in Houdini's barycentric coords
    // A - {0,0,0}; B - {1,0,0}; C - {0,1,0}

    // Check intersections between triangle 1's edges and triangle 0
    i = intersect_segment_triangle(A1, B1, C1, A0, B0, pos, parameter, uvw);
    if (i){
            append(positions, pos);
            append(uvw1, uvw);
            append(uvw0, lerp({0,0,0}, {1,0,0}, parameter));
        }

    i = intersect_segment_triangle(A1, B1, C1, B0, C0, pos, parameter, uvw);
    if (i){
            append(positions, pos);
            append(uvw1, uvw);
            append(uvw0, lerp({1,0,0}, {0,1,0}, parameter));
        }

    i = intersect_segment_triangle(A1, B1, C1, C0, A0, pos, parameter, uvw);
    if (i){
            append(positions, pos);
            append(uvw1, uvw);
            append(uvw0, lerp({0,1,0}, {0,0,0}, parameter));
        }

    // Check intersections between triangle 0's edges and triangle 1
    i = intersect_segment_triangle(A0, B0, C0, A1, B1, pos, parameter, uvw);
    if (i){
            append(positions, pos);
            append(uvw0, uvw);
            append(uvw1, lerp({0,0,0}, {1,0,0}, parameter));
        }

    i = intersect_segment_triangle(A0, B0, C0, B1, C1, pos, parameter, uvw);
    if (i){
            append(positions, pos);
            append(uvw0, uvw);
            append(uvw1, lerp({1,0,0}, {0,1,0}, parameter));
        }

    i = intersect_segment_triangle(A0, B0, C0, C1, A1, pos, parameter, uvw);
    if (i){
            append(positions, pos);
            append(uvw0, uvw);
            append(uvw1, lerp({0,1,0}, {0,0,0}, parameter));
        }

}

////////////////////////////////////////////
////// Function to check if a QUAD IS PLANAR
// https://pbr-book.org/4ed/Shapes/Bilinear_Patches

// Checks if four points are coplanar.
// Returns 1 if coplanar, 0 otherwise.
function int isCoplanar(vector p00, p01, p11, p10){
    float EPS = 1e-6f;
    // Calculate the normal of the triangle defined by p00, p01, and p10
    vector n = (normalize(cross(p10 - p00, p01 - p00)));

    // Check if the fourth point p11 is on the same plane as the triangle
    if (abs(dot(normalize(p11 - p00), n)) > EPS)
        return 0;

    return 1;
}


////////////////////////////////////////////
////// SEGMENT - PLANAR QUAD
// https://www.shadertoy.com/view/XtlBDs

// uvw - barycentric coords on the quad
// parameter - parametric position [0;1] on the segment

float cross2d( vector2 a, b ) { return a.x*b.y - a.y*b.x; } // cross product on 2d
int segment_planarquad_intersect(vector p0, p1, v0, v1, v2, v3; vector P; vector uvw; float parameter)
{
    // lets make v0 the origin

    vector ro = p0; // origin
    vector rd = p1 - p0; // direction

    vector a = v1 - v0;
    vector b = v3 - v0;
    vector c = v2 - v0;
    vector p = ro - v0;
    int lut[] = array(1,2,0,1);
    
    // intersect plane
    vector nor = cross(a,b);
    float t = -dot(p,nor) / dot(rd,nor);
    if( t<0.0 || t>1.0 ) return 0;
    
    // intersection point
    vector pos = p + t*rd;
    
    // select projection plane
    vector mor = abs(nor);
    int id = (mor.x>mor.y && mor.x>mor.z ) ? 0 : 
             (mor.y>mor.z)                 ? 1 : 
                                             2 ;

    int idu = lut[id  ];
    int idv = lut[id+1];
    
    // project to 2D
    vector2 kp = set( pos[idu], pos[idv] );
    vector2 ka = set( a[idu], a[idv] );
    vector2 kb = set( b[idu], b[idv] );
    vector2 kc = set( c[idu], c[idv] );
    
    // find barycentric coords of the quadrilateral
    vector2 kg = kc-kb-ka;

    float k0 = cross2d( kp, kb );
    float k2 = cross2d( kc-kb, ka );        // float k2 = cross2d( kg, ka );
    float k1 = cross2d( kp, kg ) - nor[id]; // float k1 = cross2d( kb, ka ) + cross2d( kp, kg );
    
    // if edges are parallel, this is a linear equation
        float u, v;
    if (abs(k2) < 0.00001){
        v = -k0/k1;
        u = cross2d( kp, ka )/k1;
        
    }else{
        // otherwise, it's a quadratic
        float w = k1*k1 - 4.0*k0*k2;
        if(w < 0.0) return 0;

        w = sqrt( w );

        float ik2 = 1.0/(2.0*k2);

                             v = (-k1 - w)*ik2; 
        if( v<0.0 || v>1.0 ) v = (-k1 + w)*ik2;
        
        u = (kp.y - ka.y * v) / (kb.y + kg.y * v);
    }
    
    if (u<0.0 || u>1.0 || v<0.0 || v>1.0) return 0;
    
    
    uvw = set(u, v, 0);
    P = v0*(1-u)*(1-v) + v1*(1-u)*v + v2*u*v + v3*u*(1-v); // houdini's function of bilinear interpolation
    parameter = t; // t - parametric position on the segment 
    
    return 1;
}



////////////////////////////////////////////
////// SEGMENT - BILINEAR PATCH (NON-PLANAR AND PLANAR QUAD)
// https://pbr-book.org/4ed/Shapes/Bilinear_Patches

// p00 - vertex 0
// p01 - vertex 1
// p11 - vertex 2
// p10 - vertex 3

int segment_bilinearpatch_intersect(vector q0, q1; vector p00, p01, p11, p10; vector P, uvw; float parameter) {
    // Checks for intersection between a segment (defined by q0 and q1) and a bilinear patch
    // defined by points p00, p01, p10, and p11.
    // Returns 1 if intersection occurs, 0 otherwise.
    // Outputs intersection point P, barycentric coordinates uvw, and intersection parameter.

    // Calculate ray parameters
    vector ray_orig = q0;
    vector ray_dir = normalize(q1 - q0);
    float tMax = length(q1 - q0);

    // Compute quadratic coefficients for distance from ray to iso-lines
    float a = dot(cross(p10 - p00, p01 - p11), ray_dir);
    float c = dot(cross(p00 - ray_orig, ray_dir), p01 - p00);
    float b = dot(cross(p10 - ray_orig, ray_dir), p11 - p10) - (a + c);

    // Solve quadratic for bilinear patch intersection
    float u1, u2;
    int roots = solvequadratic(a, b, c, u1, u2);
    if (roots == 0)
        return 0;
           
    // Find epsilon eps to ensure that candidate  is greater than zero
    // float eps = gamma(10) *
    //     (MaxComponentValue(abs(ray_orig)) + MaxComponentValue(abs(ray_dir)) +
    //      MaxComponentValue(abs(p00))   + MaxComponentValue(abs(p10))   +
    //      MaxComponentValue(abs(p01))   + MaxComponentValue(abs(p11)));
    float eps = 1e-12;
        
    // Compute v and t for the first u intersection
    float t = tMax, u, v;
    if (0 <= u1 && u1 <= 1) {
        // Precompute common terms for  and  computation
        vector uo = lerp(p00, p10, u1);
        vector ud = lerp(p01, p11, u1) - uo;
        vector deltao = uo - ray_orig;
        vector perp = cross(ray_dir, ud);
        float p2 = length2(perp);

        // Compute matrix determinants for  and  numerators
        float v1 = determinant(set(deltao.x, ray_dir.x, perp.x,
                                     deltao.y, ray_dir.y, perp.y,
                                     deltao.z, ray_dir.z, perp.z));
        float t1 = determinant(set(deltao.x, ud.x, perp.x,
                                     deltao.y, ud.y, perp.y,
                                     deltao.z, ud.z, perp.z));

        // Set u, v, and t if intersection is valid
        if (t1 > p2 * eps && 0 <= v1 && v1 <= p2) {
            u = u1;
            v = v1 / p2;
            t = t1 / p2;
        }

    }

    // Compute v and t for the second u intersection
    if (0 <= u2 && u2 <= 1 && u2 != u1) {
        vector uo = lerp(p00, p10, u2);
        vector ud = lerp(p01, p11, u2) - uo;
        vector deltao = uo - ray_orig;
        vector perp = cross(ray_dir, ud);
        float p2 = length2(perp);
        float v2 = determinant(set(deltao.x, ray_dir.x, perp.x,
                                  deltao.y, ray_dir.y, perp.y,
                                  deltao.z, ray_dir.z, perp.z));
        float t2 = determinant(set(deltao.x, ud.x, perp.x,
                                  deltao.y, ud.y, perp.y,
                                  deltao.z, ud.z, perp.z));
        t2 /= p2;
        if (0 <= v2 && v2 <= p2 && t > t2 && t2 > eps) {
           t = t2;
           u = u2;
           v = v2 / p2;
        }
    }

    // Check intersection t against tMax and possibly return intersection
    if (t >= tMax)
       return 0; // No intersection
       
    uvw = set(u, v, 0);
    P = p00*(1-u)*(1-v) + p01*(1-u)*v + p11*u*v + p10*u*(1-v);
    parameter = t / tMax;

    return 1; // Intersection found
}


