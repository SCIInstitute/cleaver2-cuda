#ifndef CAMERA_H
#define CAMERA_H

#include <Cleaver/vec3.h>

class Camera
{
public:
    virtual ~Camera(){}

    virtual Cleaver::vec3 e() const = 0;   // eye view vector
    virtual Cleaver::vec3 t() const = 0;   // tangent vector
    virtual Cleaver::vec3 u() const = 0;   // cotangent vector
    virtual Cleaver::vec3 s() const = 0;   // scale vector

    virtual void reset() = 0;
    virtual void pan(float dx, float dy) = 0;
    virtual void rotate(float theta, float phi) = 0;
    virtual void zoom(float dz) = 0;

    virtual float* viewMatrix() = 0;
};

#endif // CAMERA3D_H
