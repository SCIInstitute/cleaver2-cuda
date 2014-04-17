#ifndef PLANEFIELD_H
#define PLANEFIELD_H

#include <Cleaver/ScalarField.h>
#include <Cleaver/BoundingBox.h>


class PlaneField : public Cleaver::FloatField
{
public:
    PlaneField(double a, double b, double c, double d);
    PlaneField(const Cleaver::vec3 &n, double d);
    PlaneField(const Cleaver::vec3 &n, const Cleaver::vec3 &p);
    PlaneField(const Cleaver::vec3 &p1, const Cleaver::vec3 &p2, const Cleaver::vec3 &p3);

    virtual double valueAt(double x, double y, double z) const;
    virtual double valueAt(const Cleaver::vec3 &x) const;

    void setBounds(const Cleaver::BoundingBox &bounds);
    virtual Cleaver::BoundingBox bounds() const;

private:
    Cleaver::BoundingBox m_bounds;
    Cleaver::vec3 n;
    double d;
};

#endif // PLANEFIELD_H
