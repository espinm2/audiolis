#ifndef _GEOMETERY_UTILS_H_
#define _GEOMETERY_UTILS_H_

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <vector>
#include <cmath>
#include "ray.h"
#include "hit.h"

class ArgParser; // only for circle_points_on_plane

typedef std::vector<std::vector<double>> doubleMatrix;
/*
 * Notice:
 *   In this header you put all geometry related math that doesn't
 *   require specific knowlege of the problem domain.
 *   Only uses glm classes + ray and hit objects.
 */

// =========================================================================
// EPSILON is a necessary evil for raytracing implementations
// The appropriate value for epsilon depends on the precision of
// the floating point calculations on your hardware **AND** on the
// overall dimensions of the scene and your camera projection matrix.
#define EPSILON 0.0001
#define square(x) ((x)*(x))

// ╔═╗╦═╗╔═╗╔╦╗╔═╗╔╦╗╦ ╦╔═╗╔═╗
// ╠═╝╠╦╝║ ║ ║ ║ ║ ║ ╚╦╝╠═╝║╣ 
// ╩  ╩╚═╚═╝ ╩ ╚═╝ ╩  ╩ ╩  ╚═╝

// Parametric equation of a 3D circle
// Given two orthagonal vec a & b, a center c, radius r, and angle
glm::vec3 parametric_circle_3d(glm::vec3 a, glm::vec3 b, 
    glm::vec3 c, float r, float theta);

// Given a center n, normal n, radius r, and number of points, and empty vec
// pts we fill with points
void circle_points_on_plane( const glm::vec3 c, const glm::vec3 n, 
    const float r, const int numberPoints, std::vector<glm::vec3> &pts, 
     double offset = 0.0);

// Given a refernece point same as above
void circle_points_on_plane_refence( 
    const glm::vec3 c,  // Center of point points
    const glm::vec3 n,  // Normal of that center
    const glm::vec3 refrence,  // Refrence pointer where to start
    const float r,  //  Radius at which to set this at
    const int numberPoints,  // Number of points to generate
    std::vector<glm::vec3> &pts,  // The place we will be storing the point
    double offset = 0.0); // Offset from the refrence point

// Given points, we project them on a sphere 
// Modified pts
void cirlce_point_on_sphere( const glm::vec3 &center, const float radius, 
  std::vector<glm::vec3> &pts);

// Compute normal given 3 pts
inline glm::vec3 compute_normal(const glm::vec3 &p1, const glm::vec3 &p2, 
    const glm::vec3 &p3);

// Check if we interects a plane
bool plane_intersect(
    const Ray &r, Hit &h, 
    const glm::vec3 &a, const glm::vec3 &b, const glm::vec3 &c, 
    bool intersect_backfacing);

// Check if we interect a triangle
bool triangle_intersect( 
    const Ray &r, Hit &h, 
    const glm::vec3 &a, const glm::vec3 &b, const glm::vec3 &c, 
    bool intersect_backfacing);

inline float DistanceBetweenTwoPoints(const glm::vec3 &p1,
    const glm::vec3 &p2);

inline float AreaOfTriangle(float a, float b, float c);

inline float AreaOfTriangle(const glm::vec3 &a, 
    const glm::vec3 &b, const glm::vec3 &c);

glm::vec3 MirrorDirection(const glm::vec3 &normal, 
    const glm::vec3 &incoming);


// Untested ... what does this even do?
glm::vec3 VectorProjectPlane(const glm::vec3 & plane_normal, const glm::vec3 & v);

// Returns the point closest to points
glm::vec3 ClosestPoint(const glm::vec3 & point, const std::vector<glm::vec3>  & points);


// Returns vector of repulive forces
glm::vec3 CalcRepulsiveForces(const glm::vec3 & p, const glm::vec3 & q, 
    float r, float k);
#endif
