#ifndef SPHERESIZINGFIELD_H
#define SPHERESIZINGFIELD_H

#include <Cleaver/ScalarField.h>
#include "SphereField.h"

class SphereSizingField : public Cleaver::FloatField
{
public:
    SphereSizingField(const SphereField* field);

    virtual double valueAt(double x, double y, double z) const;
    virtual double valueAt(const Cleaver::vec3 &x) const;

    virtual Cleaver::BoundingBox bounds() const;

private:
    const SphereField *m_sphereField;
};

#endif // SPHERESIZINGFIELD_H

