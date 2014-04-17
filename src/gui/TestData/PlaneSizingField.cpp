#include "PlaneSizingField.h"
#include <algorithm>
#include <cmath>


PlaneSizingField::PlaneSizingField(const PlaneField *field)
{
    m_planeField = field;
}

double PlaneSizingField::valueAt(double x, double y, double z) const
{
    return valueAt(Cleaver::vec3(x,y,z));
}

double PlaneSizingField::valueAt(const Cleaver::vec3 &x) const
{
    double d = fabs(m_planeField->valueAt(x));

    d = std::max(0.5, d);

    return d;
}

Cleaver::BoundingBox PlaneSizingField::bounds() const
{
    return m_planeField->bounds();
}
