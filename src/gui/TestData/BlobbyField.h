#ifndef BLOBBYFIELD_H
#define BLOBBYFIELD_H

#include <Cleaver/ScalarField.h>
#include <Cleaver/BoundingBox.h>
#include <vector>

#include "SphereField.h"

class BlobbyField : public Cleaver::FloatField
{
public:
    BlobbyField(const Cleaver::vec3 &cx, float r, int nspheres, const Cleaver::BoundingBox &bounds);
    virtual double valueAt(double x, double y, double z) const;
    virtual double valueAt(const Cleaver::vec3 &x) const;

    void setBounds(const Cleaver::BoundingBox &bounds);
    virtual Cleaver::BoundingBox bounds() const;
    float meta(float r) const;
    std::vector<Cleaver::vec3> spheres;
    std::vector<float> rad;
    float a,b,c;
private:
    Cleaver::BoundingBox m_bounds;
    Cleaver::vec3 m_cx;
    float m_r;
};

#endif // SPHEREFIELD_H
