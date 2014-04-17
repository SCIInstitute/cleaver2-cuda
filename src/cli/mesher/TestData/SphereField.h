#ifndef SPHEREFIELD_H
#define SPHEREFIELD_H

#include <Cleaver/ScalarField.h>
#include <Cleaver/BoundingBox.h>

class SphereSizingField;

class SphereField : public Cleaver::FloatField
{
public:
    SphereField(const Cleaver::vec3 &cx, float r, const Cleaver::BoundingBox &bounds);

    virtual double valueAt(double x, double y, double z) const;
    virtual double valueAt(const Cleaver::vec3 &x) const;

    void setBounds(const Cleaver::BoundingBox &bounds);
    virtual Cleaver::BoundingBox bounds() const;

    friend class SphereSizingField;
    friend class SphereVaryingField;

private:
    Cleaver::BoundingBox m_bounds;
    Cleaver::vec3 m_cx;
    float m_r;
};

#endif // SPHEREFIELD_H
