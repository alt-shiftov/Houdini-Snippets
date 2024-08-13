// This is a collection of functions based on Geometric Tools. 
// I have rewritten some of these functions in VEX for use in Houdini. 
// You'll also find an example .hip file demonstrating how to use these functions.

// https://www.geometrictools.com/
// https://github.com/davideberly/GeometricTools/

/////////////////////////////
// DISTANCE POINT TO SEGMENT

function void dist_point_segment(vector point, P0, P1;
    vector out_closest[]; float out_parameter; float out_distance, out_sqrDistance){
    
    resize(out_closest, 2);

    vector direction = P1 - P0;
    vector diff = point - P1;

    float t = dot(direction, diff);
    if (t >= 0.0)
    {
        out_parameter = 1.0;
        out_closest[1] = P1;
    }
    else
    {
        diff = point - P0;
        t = dot(direction, diff);
        if (t <= 0.0)
        {
            out_parameter = 0.0;
            out_closest[1] = P0;
        }
        else
        {
            float sqrLength = dot(direction, direction);
            if (sqrLength > 0.0)
            {
                t /= sqrLength;
                out_parameter = t;
                out_closest[1] = P0 + t * direction;
            }
            else
            {
                out_parameter = 0.0;
                out_closest[1] = P0;
            }
        }
    }
    out_closest[0] = point;

    diff = out_closest[0] - out_closest[1];
    out_sqrDistance = dot(diff, diff);
    out_distance = sqrt(out_sqrDistance);

}


/////////////////////////////
// DISTANCE SEGMENT TO SEGMENT

function void dist_segment_segment(vector P0, P1, Q0, Q1; 
    vector out_closest[]; float out_parameter[]; float out_distance, out_sqrDistance){
    resize(out_closest, 2);
    resize(out_parameter, 2);

    vector P1mP0 = P1 - P0;
    vector Q1mQ0 = Q1 - Q0;
    vector P0mQ0 = P0 - Q0;
    float a = dot(P1mP0, P1mP0);
    float b = dot(P1mP0, Q1mQ0);
    float c = dot(Q1mQ0, Q1mQ0);
    float d = dot(P1mP0, P0mQ0);
    float e = dot(Q1mQ0, P0mQ0);
    float det = a * c - b * b;
    float s, t, nd, bmd, bte, ctd, bpe, ate, btd;

    if (det > 0.0){

        bte = b * e;
        ctd = c * d;
        if (bte <= ctd)  // s <= 0
        {
            s = 0.0;
            if (e <= 0.0)  // t <= 0
            {
                // region 6
                t = 0.0;
                nd = -d;
                if (nd >= a)
                {
                    s = 1.0;
                }
                else if (nd > 0.0)
                {
                    s = nd / a;
                }
                // else: s is already 0.0
            }
            else if (e < c)  // 0 < t < 1
            {
                // region 5
                t = e / c;
            }
            else  // t >= 1
            {
                // region 4
                t = 1.0;
                bmd = b - d;
                if (bmd >= a)
                {
                    s = 1.0;
                }
                else if (bmd > 0.0)
                {
                    s = bmd / a;
                }
                // else:  s is already 0.0
            }
        }
        else  // s > 0
        {
            s = bte - ctd;
            if (s >= det)  // s >= 1
            {
                // s = 1
                s = 1.0;
                bpe = b + e;
                if (bpe <= 0.0)  // t <= 0
                {
                    // region 8
                    t = 0.0;
                    nd = -d;
                    if (nd <= 0.0)
                    {
                        s = 0.0;
                    }
                    else if (nd < a)
                    {
                        s = nd / a;
                    }
                    // else: s is already one
                }
                else if (bpe < c)  // 0 < t < 1
                {
                    // region 1
                    t = bpe / c;
                }
                else  // t >= 1
                {
                    // region 2
                    t = 1.0;
                    bmd = b - d;
                    if (bmd <= 0.0)
                    {
                        s = 0.0;
                    }
                    else if (bmd < a)
                    {
                        s = bmd / a;
                    }
                    // else:  s is already one
                }
            }
            else  // 0 < s < 1
            {
                ate = a * e;
                btd = b * d;
                if (ate <= btd)  // t <= 0
                {
                    // region 7
                    t = 0.0;
                    nd = -d;
                    if (nd <= 0.0)
                    {
                        s = 0.0;
                    }
                    else if (nd >= a)
                    {
                        s = 1.0;
                    }
                    else
                    {
                        s = nd / a;
                    }
                }
                else  // t > 0
                {
                    t = ate - btd;
                    if (t >= det)  // t >= 1
                    {
                        // region 3
                        t = 1.0;
                        bmd = b - d;
                        if (bmd <= 0.0)
                        {
                            s = 0.0;
                        }
                        else if (bmd >= a)
                        {
                            s = 1.0;
                        }
                        else
                        {
                            s = bmd / a;
                        }
                    }
                    else  // 0 < t < 1
                    {
                        // region 0
                        s /= det;
                        t /= det;
                    }
                }
            }
        }
    }
    else
    {
        // The segments are parallel. The quadratic factors to
        //   R(s,t) = a*(s-(b/a)*t)^2 + 2*d*(s - (b/a)*t) + f
        // where a*c = b^2, e = b*d/a, f = |P0-Q0|^2, and b is not
        // zero. R is constant along lines of the form s-(b/a)*t = k
        // and its occurs on the line a*s - b*t + d = 0. This line
        // must intersect both the s-axis and the t-axis because 'a'
        // and 'b' are not zero. Because of parallelism, the line is
        // also represented by -b*s + c*t - e = 0.
        //
        // The code determines an edge of the domain [0,1]^2 that
        // intersects the minimum line, or if none of the edges
        // intersect, it determines the closest corner to the minimum
        // line. The conditionals are designed to test first for
        // intersection with the t-axis (s = 0) using
        // -b*s + c*t - e = 0 and then with the s-axis (t = 0) using
        // a*s - b*t + d = 0.

        // When s = 0, solve c*t - e = 0 (t = e/c).
        if (e <= 0.0)  // t <= 0
        {
            // Now solve a*s - b*t + d = 0 for t = 0 (s = -d/a).
            t = 0.0;
            nd = -d;
            if (nd <= 0.0)  // s <= 0
            {
                // region 6
                s = 0.0;
            }
            else if (nd >= a)  // s >= 1
            {
                // region 8
                s = 1.0;
            }
            else  // 0 < s < 1
            {
                // region 7
                s = nd / a;
            }
        }
        else if (e >= c)  // t >= 1
        {
            // Now solve a*s - b*t + d = 0 for t = 1 (s = (b-d)/a).
            t = 1.0;
            bmd = b - d;
            if (bmd <= 0.0)  // s <= 0
            {
                // region 4
                s = 0.0;
            }
            else if (bmd >= a)  // s >= 1
            {
                // region 2
                s = 1.0;
            }
            else  // 0 < s < 1
            {
                // region 3
                s = bmd / a;
            }
        }
        else  // 0 < t < 1
        {
            // The point (0,e/c) is on the line and domain, so we have
            // one point at which R is a minimum.
            s = 0.0;
            t = e / c;
        }
    }


    out_parameter[0] = s;
    out_parameter[1] = t;
    out_closest[0] = P0 + s * P1mP0;
    out_closest[1] = Q0 + t * Q1mQ0;
    vector diff = out_closest[0] - out_closest[1];
    out_sqrDistance = dot(diff, diff);
    out_distance = sqrt(out_sqrDistance);

}



/////////////////////////////
// DISTANCE LINE TO LINE

function void dist_line_line(vector line0_orig, line0_dir, line1_orig, line1_dir; 
    vector out_closest[]; float out_parameter[]; float out_distance, out_sqrDistance){
    resize(out_closest, 2);

    vector diff = line0_orig - line1_orig;
    float a00 = dot(line0_dir, line0_dir);
    float a01 = -dot(line0_dir, line1_dir);
    float a11 = dot(line1_dir, line1_dir);
    float b0 = dot(line0_dir, diff);
    float det = max(a00 * a11 - a01 * a01, 0.0);
    float s0, s1;
     
    if (det > 0.0){
        // The lines are not parallel.
        float b1 = -dot(line1_dir, diff);
        s0 = (a01 * b1 - a11 * b0) / det;
        s1 = (a01 * b0 - a00 * b1) / det;
    }
    else
    {
        // The lines are parallel. Select any pair of closest points.
        s0 = -b0 / a00;
        s1 = 0.0;
    }

    out_parameter[0] = s0;
    out_parameter[1] = s1;
    out_closest[0] = line0_orig + s0 * line0_dir;
    out_closest[1] = line1_orig + s1 * line1_dir;
    diff = out_closest[0] - out_closest[1];
    out_sqrDistance = dot(diff, diff);
    out_distance = sqrt(out_sqrDistance);

}


/////////////////////////////
// DISTANCE LINE TO SEGMENT

function void dist_line_segment(vector lineorig, linedir, p0, p1; 
    vector out_closest[]; float out_parameter[]; float out_distance; float out_sqrDistance){
    
    
    vector segDirection = p1 - p0;
    vector diff = lineorig - p0;
    float a00 = dot(linedir, linedir);
    float a01 = -dot(linedir, segDirection);
    float a11 = dot(segDirection, segDirection);
    float b0 = dot(linedir, diff);

    float det = max(a00 * a11 - a01 * a01, 0.0);

    float s0, s1;

    if (det > 0.0){

        // The line and segment are not parallel.
        float b1 = -dot(segDirection, diff);
        s1 = a01 * b0 - a00 * b1;

        if (s1 >= 0.0)
        {
            if (s1 <= det)
            {
                // Two interior points are closest, one on the line
                // and one on the segment.
                s0 = (a01 * b1 - a11 * b0) / det;
                s1 /= det;
            }
            else
            {
                // The endpoint Q1 of the segment and an interior
                // point of the line are closest.
                s0 = -(a01 + b0) / a00;
                s1 = 1.0;
            }
        }
        else
        {
            // The endpoint Q0 of the segment and an interior point
            // of the line are closest.
            s0 = -b0 / a00;
            s1 = 0.0;
        }
    }
    else
    {
        // The line and segment are parallel. Select the pair of
        // closest points where the closest segment point is the
        // endpoint Q0.
        s0 = -b0 / a00;
        s1 = 0.0;
    }

    // vector out_parameter;
    // vector out_closest[];
    // float out_sqrDistance;
    // float out_distance;
    
    resize(out_closest, 2);
    out_parameter[0] = s0; // длина от lineorig до точки
    out_parameter[1] = s1; // длина от p0 до точки (barycentric для сегмента)
    
    out_closest[0] = lineorig + s0 * linedir; // [0]
    out_closest[1] = p0 + s1 * segDirection; // [1]

    diff = out_closest[0] - out_closest[1];
    out_sqrDistance = dot(diff, diff);
    out_distance = sqrt(out_sqrDistance);

}


/////////////////////////////
// DISTANCE LINE TO TRIANGLE

int dist_line_triangle(vector lineorig, linedir, v0, v1, v2; 
    vector out_closest[]; float out_parameter; vector out_barycentric; float out_distance; float out_sqrDistance){
    
    
    vector E1 = v1 - v0;
    vector E2 = v2 - v0;
    vector N = cross(E1, E2);

    out_distance = 0.0;
    out_sqrDistance = 0.0;
    out_parameter = 0.0;
    out_barycentric = 0;
    out_closest = array();
    resize(out_closest, 2);
    
    float NdD = dot(N, linedir);
    if (abs(NdD) > 0.001){
        // The line and triangle are not parallel, so the line
        // intersects the plane of the triangle at a point Y.
        // Determine whether Y is contained by the triangle.
        vector PmV0 = lineorig - v0;
        float NdDiff = dot(N, PmV0);
        float tIntersect = -NdDiff / NdD;
        vector Y = lineorig + tIntersect * linedir;
        vector YmV0 = Y - v0;

        // Compute the barycentric coordinates of the intersection.
        float E1dE1 = dot(E1, E1);
        float E1dE2 = dot(E1, E2);
        float E2dE2 = dot(E2, E2);
        float E1dYmV0 = dot(E1, YmV0);
        float E2dYmV0 = dot(E2, YmV0);
        float det = E1dE1 * E2dE2 - E1dE2 * E1dE2;
        float b1 = (E2dE2 * E1dYmV0 - E1dE2 * E2dYmV0) / det;
        float b2 = (E1dE1 * E2dYmV0 - E1dE2 * E1dYmV0) / det;
        float b0 = 1.0 - b1 - b2;

        if (b0 >= 0.0 && b1 >= 0.0 && b2 >= 0.0){
            // The point Y is contained by the triangle.
            out_distance = 0.0;
            out_sqrDistance = 0.0;
            out_parameter = tIntersect;
            out_barycentric.x = b0;
            out_barycentric.y = b1;
            out_barycentric.z = b2;

            out_closest[0] = Y;
            out_closest[1] = Y;
            
            return 1; // нужно завершить раньше времени
        }
    }

    // Either (1) the line is not parallel to the triangle and the
    // point of intersection of the line and the plane of the triangle
    // is outside the triangle or (2) the line and triangle are
    // parallel. Regardless, the closest point on the triangle is on
    // an edge of the triangle. Compare the line to all three edges
    // of the triangle. To allow for arbitrary precision arithmetic,
    // the initial distance and sqrDistance are initialized to a
    // negative number rather than a floating-point maximum value.
    // Tracking the minimum requires a small amount of extra logic.


    float invalid = -1.0;
    out_distance = invalid;
    out_sqrDistance = invalid;
    
    vector triangle[] = array(v0, v1, v2);
    
 
    for (int i0 = 2, i1 = 0, i2 = 1; i1 < 3; i2 = i0, i0 = i1++){
        vector p0 = triangle[i0];
        vector p1 = triangle[i1]; 
         
        float lsOutput_sqrDistance;
        float lsOutput_distance;
        float lsOutput_parameter[];
        vector lsOutput_closest[];
        
        dist_line_segment(lineorig, linedir, p0, p1, lsOutput_closest, lsOutput_parameter,
        lsOutput_distance, lsOutput_sqrDistance);
        
        if (out_sqrDistance == invalid || lsOutput_sqrDistance < out_sqrDistance){
            out_distance = lsOutput_distance;
            out_sqrDistance = lsOutput_sqrDistance;
            out_parameter = lsOutput_parameter[0];
            out_barycentric[i0] = 1.0 - lsOutput_parameter[1];
            out_barycentric[i1] = lsOutput_parameter[1];
            out_barycentric[i2] = 0.0;
            out_closest = lsOutput_closest;
        }
    }
    return 1;
}

/////////////////////////////
// DISTANCE POINT TO TRIANGLE. Method 1

function void dist_pt_triangle(vector point, v0, v1, v2; 
    vector out_closest[]; vector out_barycentric; float out_distance; float out_sqrDistance){

    vector diff = v0 - point;
    vector edge0 = v1 - v0;
    vector edge1 = v2 - v0;

    float a00 = dot(edge0, edge0);
    float a01 = dot(edge0, edge1);
    float a11 = dot(edge1, edge1);
    float b0 = dot(diff, edge0);
    float b1 = dot(diff, edge1);
    float det = max(a00 * a11 - a01 * a01, 0.0);
    float s = a01 * b1 - a11 * b0;
    float t = a01 * b0 - a00 * b1;

    if (s + t <= det)
    {
        if (s < 0.0)
        {
            if (t < 0.0)  // region 4
            {
                if (b0 < 0.0)
                {
                    t = 0.0;
                    if (-b0 >= a00)
                    {
                        s = 1.0;
                    }
                    else
                    {
                        s = -b0 / a00;
                    }
                }
                else
                {
                    s = 0.0;
                    if (b1 >= 0.0)
                    {
                        t = 0.0;
                    }
                    else if (-b1 >= a11)
                    {
                        t = 1.0;
                    }
                    else
                    {
                        t = -b1 / a11;
                    }
                }
            }
            else  // region 3
            {
                s = 0.0;
                if (b1 >= 0.0)
                {
                    t = 0.0;
                }
                else if (-b1 >= a11)
                {
                    t = 1.0;
                }
                else
                {
                    t = -b1 / a11;
                }
            }
        }
        else if (t < 0.0)  // region 5
        {
            t = 0.0;
            if (b0 >= 0.0)
            {
                s = 0.0;
            }
            else if (-b0 >= a00)
            {
                s = 1.0;
            }
            else
            {
                s = -b0 / a00;
            }
        }
        else  // region 0
        {
            // minimum at interior point
            s /= det;
            t /= det;
        }
    }
    else
    {
        float tmp0, tmp1, numer, denom;

        if (s < 0.0)  // region 2
        {
            tmp0 = a01 + b0;
            tmp1 = a11 + b1;
            if (tmp1 > tmp0)
            {
                numer = tmp1 - tmp0;
                denom = a00 - 2.0 * a01 + a11;
                if (numer >= denom)
                {
                    s = 1.0;
                    t = 0.0;
                }
                else
                {
                    s = numer / denom;
                    t = 1.0 - s;
                }
            }
            else
            {
                s = 0.0;
                if (tmp1 <= 0.0)
                {
                    t = 1.0;
                }
                else if (b1 >= 0.0)
                {
                    t = 0.0;
                }
                else
                {
                    t = -b1 / a11;
                }
            }
        }
        else if (t < 0.0)  // region 6
        {
            tmp0 = a01 + b1;
            tmp1 = a00 + b0;
            if (tmp1 > tmp0)
            {
                numer = tmp1 - tmp0;
                denom = a00 - 2.0 * a01 + a11;
                if (numer >= denom)
                {
                    t = 1.0;
                    s = 0.0;
                }
                else
                {
                    t = numer / denom;
                    s = 1.0 - t;
                }
            }
            else
            {
                t = 0.0;
                if (tmp1 <= 0.0)
                {
                    s = 1.0;
                }
                else if (b0 >= 0.0)
                {
                    s = 0.0;
                }
                else
                {
                    s = -b0 / a00;
                }
            }
        }
        else  // region 1
        {
            numer = a11 + b1 - a01 - b0;
            if (numer <= 0.0)
            {
                s = 0.0;
                t = 1.0;
            }
            else
            {
                denom = a00 - 2.0 * a01 + a11;
                if (numer >= denom)
                {
                    s = 1.0;
                    t = 0.0;
                }
                else
                {
                    s = numer / denom;
                    t = 1.0 - s;
                }
            }
        }
    }
    resize(out_closest, 2);
    
    out_closest[0] = point; // [0]
    out_closest[1] = v0 + s * edge0 + t * edge1; // [1]

    diff = out_closest[0] - out_closest[1];
    out_sqrDistance = dot(diff, diff);
    out_distance = sqrt(out_sqrDistance);
    out_barycentric.x = 1.0 - s - t;
    out_barycentric.y = s;
    out_barycentric.z = t;

}


///////////////////////////////////////////////////////////
// DISTANCE POINT TO TRIANGLE. Method 2 - Conjugate Gradient

void GetMinEdge02(float a11, b1; float p[]){
    p[0] = 0.0;
    if (b1 >= 0.0)
    {
        p[1] = 0.0;
    }
    else if (a11 + b1 <= 0.0)
    {
        p[1] = 1.0;
    }
    else
    {
        p[1] = -b1 / a11;
    }
}

void GetMinEdge12(float a01, a11, b1, f10, f01; float p[])
{
    float h0 = a01 + b1 - f10;
    if (h0 >= 0.0)
    {
        p[1] = 0.0;
    }
    else
    {
        float h1 = a11 + b1 - f01;
        if (h1 <= 0.0)
        {
            p[1] = 1.0;
        }
        else
        {
            p[1] = h0 / (h0 - h1);
        }
    }
    p[0] = 1.0 - p[1];
}

void GetMinInterior(float p0[]; float h0; float p1[]; float h1; float p[]){
    float z = h0 / (h0 - h1);
    float omz = 1.0 - z;
    p[0] = omz * p0[0] + z * p1[0];
    p[1] = omz * p0[1] + z * p1[1];
}


function int UseConjugateGradient(vector point, v0, v1, v2; 
    vector out_closest[]; vector out_barycentric; float out_distance; float out_sqrDistance){
    
    // Compute vectors from the triangle's vertices
    vector diff = point - v0;
    vector edge0 = v1 - v0;
    vector edge1 = v2 - v0;
    
    // Compute dot products of the edge vectors with themselves and each other
    float a00 = dot(edge0, edge0);
    float a01 = dot(edge0, edge1);
    float a11 = dot(edge1, edge1);
    
    // Compute dot products of the edge vectors with the vector diff
    float b0 = -dot(diff, edge0);
    float b1 = -dot(diff, edge1);
    
    // Initialize f values, which are used for decision making in the conjugate gradient method
    float f00 = b0;
    float f10 = b0 + a00;
    float f01 = b0 + a01;
    
    // Initialize the potential solutions (p0, p1, and p)
    float p0[], p1[], p[];
    float dt1, h0, h1;

    // Compute the endpoints p0 and p1 of the segment. The segment is
    // parameterized by L(z) = (1-z)*p0 + z*p1 for z in [0,1] and the
    // directional derivative of half the quadratic on the segment is
    // H(z) = Dot(p1-p0,gradient[Q](L(z))/2), where gradient[Q]/2 =
    // (F,G). By design, F(L(z)) = 0 for cases (2), (4), (5), and
    // (6). Cases (1) and (3) can correspond to no-intersection or
    // intersection of F = 0 with the triangle.

    // Case 1: f00 >= 0, meaning the initial condition is within or on the edge of the triangle
    if (f00 >= 0.0)
    {
        // Sub-case 1.1: f01 >= 0, closest point is on the edge (0,1)
        if (f01 >= 0.0)
        {
            // (1) p0 = (0,0), p1 = (0,1), H(z) = G(L(z))
            // Find the minimum point on the edge v0-v2
            GetMinEdge02(a11, b1, p);
        }
        else
        {
            // Sub-case 1.2: Handle the interior region between p0=(0,t10) and p1=(t01,1-t01)
            
            // (2) p0 = (0,t10), p1 = (t01,1-t01),
            // H(z) = (t11 - t10)*G(L(z))
            p0[0] = 0.0;
            p0[1] = f00 / (f00 - f01);
            p1[0] = f01 / (f01 - f10);
            p1[1] = 1.0 - p1[0];
            dt1 = p1[1] - p0[1];
            h0 = dt1 * (a11 * p0[1] + b1);
            
            // Check if the minimum is on the edge or interior
            if (h0 >= 0.0)
            {
                GetMinEdge02(a11, b1, p);
            }
            else
            {
                h1 = dt1 * (a01 * p1[0] + a11 * p1[1] + b1);
                if (h1 <= 0.0)
                {
                    GetMinEdge12(a01, a11, b1, f10, f01, p);
                }
                else
                {
                    GetMinInterior(p0, h0, p1, h1, p);
                }
            }
        }
    }
    // Case 2: f01 <= 0, meaning the closest point could be on the edge v1-v2
    else if (f01 <= 0.0)
    {
        // Sub-case 2.1: f10 <= 0, closest point is on the edge (1,0)
        if (f10 <= 0.0)
        {
            // (3) p0 = (1,0), p1 = (0,1), H(z) = G(L(z)) - F(L(z))
            GetMinEdge12(a01, a11, b1, f10, f01, p);
        }
        
        // Sub-case 2.2: Handle the interior region between p0=(t00,0) and p1=(t01,1-t01)
        else
        {
            // (4) p0 = (t00,0), p1 = (t01,1-t01), H(z) = t11*G(L(z))
            p0[0] = f00 / (f00 - f10);
            p0[1] = 0.0;
            p1[0] = f01 / (f01 - f10);
            p1[1] = 1.0 - p1[0];
            h0 = p1[1] * (a01 * p0[0] + b1);
            
            // Check if the minimum is on the edge or interior
            if (h0 >= 0.0)
            {
                p = p0;  // // Use p0 as the closest point
            }
            else
            {
                h1 = p1[1] * (a01 * p1[0] + a11 * p1[1] + b1);
                if (h1 <= 0.0)
                {
                    GetMinEdge12(a01, a11, b1, f10, f01, p);
                }
                else
                {
                    GetMinInterior(p0, h0, p1, h1, p);
                }
            }
        }
    }
    // Case 3: f10 <= 0, meaning the closest point could be on the edge v0-v1
    else if (f10 <= 0.0)
    {
        // Sub-case 3.1: Handle the interior region between p0=(0,t10) and p1=(t01,1-t01)
        // (5) p0 = (0,t10), p1 = (t01,1-t01),
        // H(z) = (t11 - t10)*G(L(z))
        p0[0] = 0.0;
        p0[1] = f00 / (f00 - f01);
        p1[0] = f01 / (f01 - f10);
        p1[1] = 1.0 - p1[0];
        dt1 = p1[1] - p0[1];
        h0 = dt1 * (a11 * p0[1] + b1);
        
        // Check if the minimum is on the edge or interior
        if (h0 >= 0.0)
        {
            GetMinEdge02(a11, b1, p);
        }
        else
        {
            h1 = dt1 * (a01 * p1[0] + a11 * p1[1] + b1);
            if (h1 <= 0.0)
            {
                GetMinEdge12(a01, a11, b1, f10, f01, p);
            }
            else
            {
                GetMinInterior(p0, h0, p1, h1, p);
            }
        }
    }
    // Case 4: Handle the general case where the closest point is inside the triangle
    else
    {
        // Sub-case 4.1: Handle the interior region between p0=(t00,0) and p1=(0,t11)
        // (6) p0 = (t00,0), p1 = (0,t11), H(z) = t11*G(L(z))
        p0[0] = f00 / (f00 - f10);
        p0[1] = 0.0;
        p1[0] = 0.0;
        p1[1] = f00 / (f00 - f01);
        h0 = p1[1] * (a01 * p0[0] + b1);
        
        // Check if the minimum is on the edge or interior
        if (h0 >= 0.0)
        {
            p = p0;  // Use p0 as the closest point
        }
        else
        {
            h1 = p1[1] * (a11 * p1[1] + b1);
            if (h1 <= 0.0)
            {
                GetMinEdge02(a11, b1, p);
            }
            else
            {
                GetMinInterior(p0, h0, p1, h1, p);
            }
        }
    }

    append(out_closest, point); // [0]
    append(out_closest, v0 + p[0] * edge0 + p[1] * edge1); // [1]

    diff = out_closest[0] - out_closest[1];
    out_sqrDistance = dot(diff, diff);
    out_distance = sqrt(out_sqrDistance);
    
    out_barycentric.x = 1.0 - p[0] - p[1];
    out_barycentric.y = p[0];
    out_barycentric.z = p[1];

    return 1;

}



///////////////////////////////
// SEGMENT TO TRIANGLE
///////////////////////////////

function void dist_segment_triangle(vector p0, p1, v0, v1, v2; 
    vector out_closest[]; float out_parameter; vector out_barycentric; float out_distance; float out_sqrDistance){
    
    resize(out_closest, 2);
    
    vector segDirection = p1 - p0;
    vector lineorig = p0;
    vector linedir = (segDirection);

    float ltOutput_parameter; vector ltOutput_closest[];
    vector ltOutput_barycentric; float ltOutput_distance; float ltOutput_sqrDistance;
    
    vector ptOutput_closest[]; vector ptOutput_barycentric; 
    float ptOutput_distance; float ptOutput_sqrDistance;    
    
    dist_line_triangle(lineorig, linedir, v0, v1, v2, ltOutput_closest, ltOutput_parameter, ltOutput_barycentric, ltOutput_distance, ltOutput_sqrDistance);
    
    if (ltOutput_parameter >= 0.0){
        
        if (ltOutput_parameter <= 1.0){
            out_closest = ltOutput_closest;
            out_parameter = ltOutput_parameter;
            out_barycentric = ltOutput_barycentric;
            out_distance = ltOutput_distance;
            out_sqrDistance = ltOutput_sqrDistance;
        }else{
        
            dist_pt_triangle(p1, v0, v1, v2, ptOutput_closest, ptOutput_barycentric, ptOutput_distance, ptOutput_sqrDistance);
            
            out_distance = ptOutput_distance;
            out_sqrDistance = ptOutput_sqrDistance;
            out_parameter = 1.0;
            out_barycentric = ptOutput_barycentric;
            out_closest[0] = p1;
            out_closest[1] = ptOutput_closest[1];
        }
    }else{

        dist_pt_triangle(p0, v0, v1, v2, ptOutput_closest, ptOutput_barycentric, ptOutput_distance, ptOutput_sqrDistance);
        
        out_distance = ptOutput_distance;
        out_sqrDistance = ptOutput_sqrDistance;
        out_parameter = 0.0;
        out_barycentric = ptOutput_barycentric;
        
        out_closest[0] = p0;
        out_closest[1] = ptOutput_closest[1];
    }

}

/////////////////////////
// TRIANGLE TO TRIANGLE
/////////////////////////

// не выдает все точки пересечения. Выдает по одной ближайшей точке на каждом треугольнике

// v0, v1, v2 - triangle 1
// w0, w1, w2 - triangle 2
function void dist_triangle_triangle(vector v0, v1, v2, w0, w1, w2; 
    vector out_closest[]; float out_parameter; vector out_barycentric0, out_barycentric1; float out_distance, out_sqrDistance){
    
    float invalid = -1.0;
    out_distance = invalid;
    out_sqrDistance = invalid;
    resize(out_closest, 2);
    
    vector triangle0[] = array(v0, v1, v2);
    vector triangle1[] = array(w0, w1, w2);
    
    float lsOutput_sqrDistance;
    float lsOutput_distance;
    vector lsOutput_barycentric;
    float lsOutput_parameter;
    vector lsOutput_closest[];
    
    // Compare edges of triangle0 to the interior of triangle1.
    for (int i0 = 2, i1 = 0, i2 = 1; i1 < 3; i2 = i0, i0 = i1++){
        vector p0 = triangle0[i0];
        vector p1 = triangle0[i1];
        
        // float lsOutput_sqrDistance;
        // float lsOutput_distance;
        // vector lsOutput_barycentric;
        // float lsOutput_parameter;
        // vector lsOutput_closest[];

        dist_segment_triangle(p0, p1, w0, w1, w2, lsOutput_closest, lsOutput_parameter, lsOutput_barycentric, lsOutput_distance, lsOutput_sqrDistance);
        
        if (out_sqrDistance == invalid || lsOutput_sqrDistance < out_sqrDistance){
            out_distance = lsOutput_distance;
            out_sqrDistance = lsOutput_sqrDistance;

            out_barycentric0[i0] = 1.0 - lsOutput_parameter;
            out_barycentric0[i1] = lsOutput_parameter;
            out_barycentric0[i2] = 0.0;
            
            out_barycentric1 = lsOutput_barycentric;
            out_closest = lsOutput_closest;
        }
    }

    // Compare edges of triangle1 to the interior of triangle0.
    for (int i0 = 2, i1 = 0, i2 = 1; i1 < 3; i2 = i0, i0 = i1++){
        vector p0 = triangle1[i0];
        vector p1 = triangle1[i1];
        
        dist_segment_triangle(p0, p1, v0, v1, v2, lsOutput_closest, lsOutput_parameter, lsOutput_barycentric, lsOutput_distance, lsOutput_sqrDistance);
        
        if (out_sqrDistance == invalid || lsOutput_sqrDistance < out_sqrDistance){
            out_distance = lsOutput_distance;
            out_sqrDistance = lsOutput_sqrDistance;
            out_barycentric0 = lsOutput_barycentric;

            out_barycentric1[i0] = 1.0 - lsOutput_parameter;
            out_barycentric1[i1] = lsOutput_parameter;
            out_barycentric1[i2] = 0.0;
            
            out_closest[0] = lsOutput_closest[1];
            out_closest[1] = lsOutput_closest[0];
        }
    }
    
}
