#include "TargetCamera.h"
#include <cmath>

using namespace Cleaver;

TargetCamera::TargetCamera()
{
    reset();
}

Cleaver::vec3 TargetCamera::e() const
{
    return m_e;
}

Cleaver::vec3 TargetCamera::t() const
{
    return m_t;
}

Cleaver::vec3 TargetCamera::u() const
{
    return m_u;
}

Cleaver::vec3 TargetCamera::s() const
{
    return Cleaver::vec3(m_scale, m_scale, m_scale);
}

float* TargetCamera::viewMatrix()
{
    return m_viewMatrix;
}

void TargetCamera::computeViewMatrix()
{

    Cleaver::vec3 forward = Cleaver::normalize(m_t - m_e);
    Cleaver::vec3    side = Cleaver::normalize(cross(forward, m_u));
    Cleaver::vec3      up = Cleaver::normalize(cross(side, forward));

    //------------------
    m_viewMatrix[0]  = side.x;
    m_viewMatrix[4]  = side.y;
    m_viewMatrix[8]  = side.z;
    m_viewMatrix[12] = 0;
    //------------------
    m_viewMatrix[1]  = up.x;
    m_viewMatrix[5]  = up.y;
    m_viewMatrix[9]  = up.z;
    m_viewMatrix[13] = 0;
    //------------------
    m_viewMatrix[2]  = -forward.x;
    m_viewMatrix[6]  = -forward.y;
    m_viewMatrix[10] = -forward.z;
    m_viewMatrix[14] = 0;
    //------------------

    m_viewMatrix[12] = -1*dot(    side, m_e);
    m_viewMatrix[13] = -1*dot(      up, m_e);
    m_viewMatrix[14] =    dot( forward, m_e);
    //------------------
    m_viewMatrix[3] = 0;
    m_viewMatrix[7] = 0;
    m_viewMatrix[11] = 0;
    m_viewMatrix[15] = 1.0;
}

void TargetCamera::reset()
{
    float view_distance = 1.5*Cleaver::length(m_targetBounds.center() - m_targetBounds.origin);
    m_t   = m_targetBounds.center();
    m_e.x = m_t.x + view_distance;
    m_e.y = m_t.y + view_distance;
    m_e.z = m_t.z +-view_distance;

    m_scale = length(m_t - m_e);

    m_u.x = 0;
    m_u.y = 1;
    m_u.z = 0;

    m_viewDir = normalize(m_t - m_e);
    m_r = cross(m_viewDir, vec3(0,1,0));
    m_u = cross(m_r, m_viewDir);

    computeViewMatrix();
}

void TargetCamera::setView(const vec3 &eye, const vec3 &target)
{
    m_e = eye;
    m_t = target;

    m_viewDir = normalize(m_t - m_e);
    m_r = cross(m_viewDir, vec3(0,1,0));
    m_u = cross(m_r, m_viewDir);

    m_scale = length(m_t - m_e);

    computeViewMatrix();
}

void TargetCamera::zoom(float dz)
{
    if(dz > 0)
        m_scale /= 1.05;
    else
        m_scale *= 1.05;

    // don't allow inversion
    if(m_scale == 0)
        m_scale += 0.0001;

    m_e = m_t - m_scale*m_viewDir;

    computeViewMatrix();
}

void TargetCamera::pan(float dx, float dy)
{
    m_e += 0.05*(dx*m_r + dy*m_u);
    m_t += 0.05*(dx*m_r + dy*m_u);
    computeViewMatrix();
}

void TargetCamera::rotate(float theta, float phi)
{
    theta *= 0.002f;
    phi   *= 0.002f;
    vec3 origin = m_e - m_t;
    m_e.z = origin.z*cos(theta) - origin.x*sin(theta);
    m_e.x = origin.z*sin(theta) + origin.x*cos(theta);
    m_e.y = origin.y;
    m_e += m_t;

    this->setView(m_e, m_t);

    phi += asin((m_e.y - m_t.y) / m_scale);
    float h = sin(phi)*m_scale;
    float d = cos(phi)*m_scale;

    vec3 flatView(m_viewDir.x, 0, m_viewDir.z);
    flatView = normalize(flatView);

    m_e.x = m_t.x - d*flatView.x;
    m_e.y = m_t.y + h;
    m_e.z = m_t.z - d*flatView.z;

    this->setView(m_e, m_t);

    computeViewMatrix();
}

void TargetCamera::setTargetBounds(const Cleaver::BoundingBox &bounds)
{
    m_targetBounds = bounds;
}
