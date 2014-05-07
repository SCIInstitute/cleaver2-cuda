#ifndef DATAMANAGER_H
#define DATAMANAGER_H

#include <Cleaver/TetMesh.h>
#include <Cleaver/Volume.h>
#include <QObject>
#include <stdlib.h>
#include <array>

template<typename T>
std::ostream& binary_write(std::ostream& stream, const T& value){
    return stream.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

class DataManager : public QObject
{
  Q_OBJECT

 public:
  DataManager();

  void addMesh(Cleaver::TetMesh *mesh);
  void removeMesh(Cleaver::TetMesh *mesh);

  void addField(Cleaver::AbstractScalarField *field);
  void removeField(Cleaver::AbstractScalarField *field);
  void clearFields();

  void addVolume(Cleaver::Volume *volume);
  void addTansitionMesh(std::vector<std::array<float,3>> &verts,
                        std::vector<std::array<size_t,4>> &faces,
                        std::array<size_t,3> dims);
  void outputTansitionMesh(std::string name);
  void removeVolume(Cleaver::Volume *volume);

  std::vector<ulong> getSelection();
  void setSelection(ulong);
  void addSelection(ulong);
  void toggleSetSelection(ulong);
  void toggleAddSelection(ulong);
  void clearSelection();

  void update(){ emit dataChanged(); }

  std::vector<Cleaver::AbstractScalarField*>  fields() const { return m_fields; }
  std::vector<Cleaver::Volume*>       volumes() const { return m_volumes; }
  std::vector<Cleaver::TetMesh*>      meshes() const { return m_meshes; }

  signals:

  void dataAdded();
  void dataRemoved();
  void dataChanged();

  void meshAdded();
  void meshRemoved();

  void fieldAdded();
  void fieldRemoved();

  void volumeAdded();
  void volumeRemoved();

  void meshListChanged();
  void fieldListChanged();
  void volumeListChanged();

  void selectionChanged();

 private:

  std::vector<ulong>                          m_selection;
  std::vector<Cleaver::TetMesh*>              m_meshes;
  std::vector<Cleaver::Volume*>               m_volumes;
  std::vector<Cleaver::AbstractScalarField*>  m_fields;
  std::vector<std::array<float,3>>            verts_;
  std::vector<std::array<size_t,4>>           faces_;
  std::array<size_t, 3>                       dims_;
};

#endif // DATAMANAGER_H
