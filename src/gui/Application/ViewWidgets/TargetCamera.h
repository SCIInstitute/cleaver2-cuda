#ifndef TARGETCAMERA_H
#define TARGETCAMERA_H

#include "Camera.h"
#include <Cleaver/vec3.h>
#include <Cleaver/BoundingBox.h>

class TargetCamera : public Camera
{
public:
    TargetCamera();

    virtual Cleaver::vec3 e() const;   // eye location
    virtual Cleaver::vec3 t() const;   // target location
    virtual Cleaver::vec3 u() const;   // up vector
    virtual Cleaver::vec3 s() const;   // scale vector

    virtual float* viewMatrix(); // return viewMatrix in Col-Major

    virtual void reset();
    virtual void pan(float dx, float dy);
    virtual void rotate(float theta, float phi);
    virtual void zoom(float dz);

    void setView(const Cleaver::vec3 &eye, const Cleaver::vec3 &target);
    void setTargetBounds(const Cleaver::BoundingBox &bounds);

private:

    void computeViewMatrix();

    Cleaver::vec3 m_e;       // eye location
    Cleaver::vec3 m_t;       // target location
    Cleaver::vec3 m_u;       // up direction
    Cleaver::vec3 m_r;       // right direction
    Cleaver::vec3 m_viewDir; // view direction
    Cleaver::BoundingBox m_targetBounds; // target bounds

    float m_scale;
    float m_viewMatrix[16];
};

#endif // TARGETCAMERA_H
