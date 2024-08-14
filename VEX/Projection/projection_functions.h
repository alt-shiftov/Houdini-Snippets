// Projection simply involves finding the minimal distance.

// This can be achieved using Near Distance Functions (Geometric Tools Functions).
// https://github.com/alt-shiftov/Houdini-Snippets/tree/main/VEX/NearestDistance
// https://github.com/davideberly/GeometricTools

// However, I found another implementations that might also be useful.


/////////////////////////////////////
// CHECK IF POINT IS INSIDE THE TRIANGLE

int point_in_triangle(vector P, A, B, C; vector uvw)
{
    float EPS = 1e-6;
    uvw = set(0,0,0);
    
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
    if (abs(denom) < EPS) {
        return 0; // Degenerate triangle
    }

    float u = (dot11 * dot02 - dot01 * dot12) / denom;
    if (u < 0.0 || u > 1.0) {
        return 0; // Point outside triangle
    }

    float v = (dot00 * dot12 - dot01 * dot02) / denom;
    if (v < 0.0 || u + v > 1.0) {
        return 0; // Point outside triangle
    }
    uvw = set(u, v, 0);
    
    return 1; // Point inside triangle
}


////////////////////////////
// PROJECT POINT ON A PLANE

vector project_point_on_plane(vector pt, plane_n, plane_pt)
{
    // Calculate the distance from point to the plane
    float distance = dot(pt - plane_pt, plane_n);

    // Project point A onto the plane
    vector projectedPoint = pt - distance * plane_n;

    return projectedPoint;
}


////////////////////////////
// PROJECT POINT ON A SEGMENT (version 1)

vector project_point_on_segment(vector A, P0, P1; float parameter)
{
    // Vector from P0 to P1
    vector seg = P1 - P0;
    
    // Vector from P0 to A
    vector AP = A - P0;
    
    // Projection length along the segment
    float t = dot(AP, seg) / dot(seg, seg);
    
    // to project a point onto a line, simply remove the clamping parameter.
    // Clamp t to [0, 1] to ensure projection is on the segment
    t = clamp(t, 0.0, 1.0);
    
    // Projected point
    vector B = P0 + t * seg;
    parameter = t;
    
    return B;
}

// PROJECT POINT ON A SEGMENT (version 2)
// https://paulbourke.net/geometry/pointlineplane/source.c

int project_point_on_segment( vector point, start, end, intersection; float parameter )
{
    float linemag = distance(end, start);
 
    float u = ( ( ( point.x - start.x ) * ( end.x - start.x ) ) +
        ( ( point.y - start.y ) * ( end.y - start.y ) ) +
        ( ( point.z - start.z ) * ( end.z - start.z ) ) ) /
        ( linemag * linemag );

    // to project a point onto a line, simply remove the clamping parameter.
    u = clamp(u, 0.0f, 1.0f);
 
    intersection.x = start.x + u * ( end.x - start.x );
    intersection.y = start.y + u * ( end.y - start.y );
    intersection.z = start.z + u * ( end.z - start.z );
    parameter = u;
    
    return 1;
}


/////////////////////////////////
// PROJECT POINT ON A TRIANGLE

// Returns the position of the projected point.
// Method may be less accurate near the corners of the triangle.
// For better results, consider using dist_pt_triangle() from neardistance.h.

function vector project_point_on_triangle(vector a, b, c, p; vector uvw) {

    vector v0 = b - a; 
    vector v1 = c - a;
    vector v2 = p - a;
    
    float d00 = dot(v0, v0);
    float d01 = dot(v0, v1);
    float d11 = dot(v1, v1);
    float d20 = dot(v2, v0);
    float d21 = dot(v2, v1);
    
    float denom = d00 * d11 - d01 * d01;
    
    float u = (d11 * d20 - d01 * d21) / denom;
    float v = (d00 * d21 - d01 * d20) / denom;
    float w = 1.0f - u - v;
    
    vector p0 = b;
    vector p1 = c;
    vector p2 = a;
    
    if ( u < 0){
            float t = dot(p-p1,p2-p1)/dot(p2-p1,p2-p1);
            t = clamp( t, 0, 1 );
            uvw = set( 0.0f, 1.0f-t, t );
            
    }else if ( v < 0 ){
            float t = dot(p-p2,p0-p2)/dot(p0-p2,p0-p2);
            t = clamp( t, 0, 1 );
            uvw = set( t, 0.0f, 1.0f-t );
    }else if ( w < 0 ){
            float t = dot(p-p0,p1-p0)/dot(p1-p0,p1-p0);
            t = clamp( t, 0, 1 );
            uvw = set( 1.0f-t, t, 0.0f );
    }else{
            uvw = set( u, v, 0 );
    }
    
    u = uvw.x;
    v = uvw.y;
    
    return a*(1-u-v) + b*u + c*v;
}


//////////////////////////////////////
// Function to check if a QUAD IS PLANAR
// https://pbr-book.org/4ed/Shapes/Bilinear_Patches

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



///////////////////////////////////////
// PROJECT POINT ON PLANAR QUAD

// This is a very clunky method.

// projectPointOnQuad() function determines the projected position and UVW coordinates.
// First, it orients the quad to 2D space using orientquad() (excluding the Y coordinate).
// Second, it finds the closest point on the quad with closestPointOnQuad().
// Finally, it computes the UVW Barycentric coordinates for the closest position using invBilinear().

// invBilinear() - by https://iquilezles.org/articles/ibilinear/

// function to orient quad and point position to XY space
function vector4 orientquad(vector p, p00, p01, p11, p10){
    // finding quat for orienting quad to XY
    vector n = (normalize(cross(p10 - p00, p01 - p00)));
    vector4 quat = dihedral(n, set(0,0,1));
    
    // orient all quad vertices and source point
    p00 = qrotate(quat, p00);
    p01 = qrotate(quat, p01);
    p11 = qrotate(quat, p11);
    p10 = qrotate(quat, p10);
    
    p = qrotate(quat, p);
    
    return quat;
}


// Function to check if point is in triangle using bary coords.
int pointInTriangle(vector p, a, b, c) {
    vector v0 = c - a;
    vector v1 = b - a;
    vector v2 = p - a;

    float dot00 = dot(v0, v0);
    float dot01 = dot(v0, v1);
    float dot02 = dot(v0, v2);
    float dot11 = dot(v1, v1);
    float dot12 = dot(v1, v2);

    float invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
    float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    return (u >= 0) && (v >= 0) && (u + v <= 1);
}

// Function for finding the closest point on the quad's edges
vector closestPointOnLine(vector p, a, b) {
    vector ab = b - a;
    float t = dot(p - a, ab) / dot(ab, ab);
    t = clamp(t, 0.0, 1.0);  // Clamp to line segment
    return a + t * ab;
}


function vector closestPointOnQuad(vector initpos, q0, q1, q2, q3) {
    vector n = (normalize(cross(q2 - q0, q1 - q0)));
    
    // Calculate the distance from point A to the plane
    float dist = dot(initpos - q0, n);

    // Project point A onto the plane
    vector p = initpos - dist * n;
    
    // Barycentric coordinates method to check if the point is inside the quad
    vector2 p2d = set(p.x, p.y);
    vector2 q02d = set(q0.x, q0.y);
    vector2 q12d = set(q1.x, q1.y);
    vector2 q22d = set(q2.x, q2.y);
    vector2 q32d = set(q3.x, q3.y);

    // Check if point is inside either of the two triangles forming the quad
    if (pointInTriangle(p2d, q02d, q12d, q22d) || pointInTriangle(p2d, q02d, q22d, q32d)) {
        return p;
    }
    // If the point is outside, find the closest point on quad's edge
    vector cp0 = closestPointOnLine(p, q0, q1);
    vector cp1 = closestPointOnLine(p, q1, q2);
    vector cp2 = closestPointOnLine(p, q2, q3);
    vector cp3 = closestPointOnLine(p, q3, q0);

    // Find the closest of the four points
    float dist0 = distance2(p, cp0);
    float dist1 = distance2(p, cp1);
    float dist2 = distance2(p, cp2);
    float dist3 = distance2(p, cp3);

    if (dist0 < dist1 && dist0 < dist2 && dist0 < dist3) {
        return cp0;
    } else if (dist1 < dist2 && dist1 < dist3) {
        return cp1;
    } else if (dist2 < dist3) {
        return cp2;
    } else {
        return cp3;
    }
}

// by IQ: https://iquilezles.org/articles/ibilinear/
float cross2d(vector a, b ) { return a.x*b.y - a.y*b.x; }
vector invBilinear(vector p, a, b, c, d )
{
    float EPS = 1e-6;
    vector res = set(-1.0, -1, -1);

    vector e = b-a;
    vector f = d-a;
    vector g = a-b+c-d;
    vector h = p-a;
        
    float k2 = cross2d( g, f );
    float k1 = cross2d( e, f ) + cross2d( h, g );
    float k0 = cross2d( h, e );
    
    // if edges are parallel, this is a linear equation
    if( abs(k2)<EPS )
    {
        res = set( (h.x*k1+f.x*k0)/(e.x*k1-g.x*k0), -k0/k1 );
    }
    // otherwise, it's a quadratic
    else
    {
        float w = k1*k1 - 4.0*k0*k2;
        // if( w<0.0 ) return set(-1.0, -1, -1);
        if (w<0.0) w = 0.0;
        
        w = sqrt( w );
        printf("%g\n", w);
        
        float ik2 = 0.5/k2;
        float v = (-k1 - w)*ik2;
        float u = (h.x - f.x*v)/(e.x + g.x*v);
        
        // if( u<0.0 || u>1.0 || v<0.0 || v>1.0 )
        // {
        //    v = (-k1 + w)*ik2;
        //    u = (h.x - f.x*v)/(e.x + g.x*v);
        // }
        res = set( v, u, 0);
    }
    
    return res;
}


function vector projectPointOnQuad(vector p, p00, p01, p11, p10, uvw){
    vector p_2d = p;
    vector p00_2d = p00;
    vector p01_2d = p01;
    vector p11_2d = p11;
    vector p10_2d = p10;
    
    // First, orient the quad  (removing Z axis)
    vector4 quat = orientquad(p_2d, p00_2d, p01_2d, p11_2d, p10_2d);
    
    // Second, project the point onto the quad (if necessary, clamping to the borders)
    vector projpos2d = closestPointOnQuad(p_2d, p00_2d, p01_2d, p11_2d, p10_2d);
    
    // Finally, calculate the UVW coordinates for the projected position.
    uvw = invBilinear(projpos2d, p00_2d, p01_2d, p11_2d, p10_2d); 
    
    // Aaand reorient the projected position back to XYZ coordinates.
    return qrotate(qinvert(quat), projpos2d);
}




///////////////////////////////////////
// PROJECT POINT ON BILINEAR PATCH (planar and non-planar quad)

// This method uses iterative coordinate-wise minimization to find the closest bary coords in a quad.
// A step value of 7-8 yields good results.
// For planar quads, the results are similar to those produced by the xyzdist() function.
// For non-planar quads, the results may differ from those of xyzdist().

// By IQ - https://www.shadertoy.com/view/3tjczm

float dot2( vector v ) { return dot(v,v); }
vector closestPointToBilinear( vector p, p0, p1, p2, p3, uvw; int step){    
    vector A = p1-p0;
    vector B = p3-p0;
    vector C = p2-p3-p1+p0;
    vector D = p0-p;

    // initial guess
    uvw = set(0.0,0.0);
    
    float d = dot2(p-p0);
    int nmax = (int)pow(2, step);
    
    for( int i=0; i<nmax; i++ ){
        float u = float(i)/((nmax-1)*1.0);
        vector ba = lerp( B,p2-p1,u);
        vector pa = lerp(-D,p -p1,u);
        float v = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
        float t = dot2(pa-ba*v);
        if( t<d ){ 
            d=t; 
            uvw = set(v,u,0);
        }
    }
    float u = uvw.x;
    float v = uvw.y;
    
    return p0*(1-u)*(1-v) + p1*(1-u)*v + p2*u*v + p3*u*(1-v);
}

/////////////////////////////////////////




