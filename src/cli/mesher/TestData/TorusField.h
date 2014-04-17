#ifndef TORUSFIELD_H
#define TORUSFIELD_H

#include <Cleaver/ScalarField.h>
#include <Cleaver/BoundingBox.h>
#include <vector>

class TorusField : public Cleaver::FloatField
{
public:
    TorusField(const Cleaver::vec3 &cx, float ur, float vr, const Cleaver::BoundingBox &bounds);

    virtual double valueAt(double x, double y, double z) const;
    virtual double valueAt(const Cleaver::vec3 &x) const;

    void setBounds(const Cleaver::BoundingBox &bounds);
    virtual Cleaver::BoundingBox bounds() const;

    std::vector<Cleaver::vec3>  tensorAt(const Cleaver::vec3 &x) const;


private:
    Cleaver::BoundingBox m_bounds;
    Cleaver::vec3 m_cx;             // center
    float m_ur;                     // minor radius
    float m_vr;                     // minor radius
};

#endif // TORUSFIELD_H
