#ifndef PLANESIZINGFIELD_H
#define PLANESIZINGFIELD_H

#include <Cleaver/ScalarField.h>
#include "PlaneField.h"

class PlaneSizingField : public Cleaver::FloatField
{
public:
    PlaneSizingField(const PlaneField *field);

    virtual double valueAt(double x, double y, double z) const;
    virtual double valueAt(const Cleaver::vec3 &x) const;

    virtual Cleaver::BoundingBox bounds() const;

private:
    const PlaneField *m_planeField;
};
#endif // PLANESIZINGFIELD_H
