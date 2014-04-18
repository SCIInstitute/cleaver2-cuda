#ifndef MESHWINDOW_H
#define MESHWINDOW_H

// enable vertex buffer includes
#define GL_GLEXT_PROTOTYPES

#if defined(WIN32)
#include <GL/glew.h>
#endif

#include <QGLWidget>
#include <vector>
#include <array>
#include <Cleaver/CleaverMesher.h>
#include <Cleaver/TetMesh.h>
#include <Cleaver/Volume.h>
#include "Camera.h"

enum StarMode { NoStar, VertexStar, EdgeStar, FaceStar };
enum CameraType { Target, Trackball };


class MeshWindow : public QGLWidget
{
  Q_OBJECT
 public:
  explicit MeshWindow(QObject *parent = 0);
  ~MeshWindow();

  bool connectToTimer(QTimer *timer);
  void setDefaultOptions();
  void resetView();
  void updateMesh();

  //-- mutators --
  void setMesh(Cleaver::TetMesh *mesh);
  void setMesh(std::vector<std::array<float,3>> &verts,
               std::vector<std::array<size_t,4>> &faces);
  void setVolume(Cleaver::Volume *volume);
  void setAxisVisible(bool value){ m_bShowAxis = value; }
  void setBBoxVisible(bool value){ m_bShowBBox = value; }
  void setParticlesVisible(bool value){ m_bShowParticles = value; }
  void setFacesVisible(bool value){ m_bShowFaces = value; }
  void setEdgesVisible(bool value){ m_bShowEdges = value; }
  void setCutsVisible(bool value){ m_bShowCuts = value; }
  void setSyncedClipping(bool value){ m_bSyncedClipping = value;  }
  void setSurfacesOnly(bool value){ m_bSurfacesOnly = value; }
  void setColorByQuality(bool value){ m_bColorByQuality = value; }
  void setClippingPlaneVisible(bool value){ m_bShowClippingPlane = value; }
  void setClipping(bool value){ m_bClipping = value; update_vbos(); }
  void setClippingPlane(float plane[4]){ memcpy(m_4fvClippingPlane, plane, 4*sizeof(float)); if(!m_bShowClippingPlane || m_bSyncedClipping) update_vbos();}
  void setCamera(Camera *camera){ m_camera = camera; }
  void setCameraType(CameraType type){ m_cameraType = type; initializeCamera(); }
  void setMaterialFaceLock(int m, bool value){ m_bMaterialFaceLock[m] = value; }
  void setMaterialCellLock(int m, bool value){ m_bMaterialCellLock[m] = value; }


  void setParticles(std::vector<Cleaver::vec3> *p){ particles = p; }
  //void setClippingPlane(float plane[4]){ memcpy(m_4fvClippingPlane, plane, 4*sizeof(float)); update_vbos();}
  //void setClipping(bool value){ m_bClipping = value; update_vbos(); }

  //-- acccessors
  bool axisVisible(){ return m_bShowAxis; }
  bool bboxVisible(){ return m_bShowBBox; }
  bool facesVisible(){ return m_bShowFaces; }
  bool edgesVisible(){ return m_bShowEdges; }
  bool cutsVisible(){ return m_bShowCuts; }
  size_t max_mat() { return max_mat_; }

  bool getMaterialFaceLock(int m) const { return m_bMaterialFaceLock[m]; }
  bool getMaterialCellLock(int m) const { return m_bMaterialCellLock[m]; }

  Cleaver::CleaverMesher* mesher(){ return m_mesher; }
  Cleaver::TetMesh* mesh(){ return m_mesh; }
  Cleaver::Volume* volume(){ return m_volume; }
  Cleaver::BoundingBox dataBounds(){ return m_dataBounds; }
  bool validMesh() { return verts_.size() > 0 && faces_.size() > 0; }


  signals:

  public slots:

  private:

  Cleaver::CleaverMesher *m_mesher;
  Cleaver::TetMesh *m_mesh;
  Cleaver::Volume  *m_volume;
  std::vector<std::array<float,3>>            verts_;
  std::vector<std::array<size_t,4>>           faces_;
  size_t max_mat_ = 0;
  Cleaver::BoundingBox m_dataBounds;
  Camera *m_camera;
  CameraType m_cameraType;


  StarMode m_starmode;
  int m_currentVertex;
  int m_currentEdge;
  int m_currentFace;


  int m_prev_x, m_prev_y;

  // render options
  bool m_bShowAxis;
  bool m_bShowBBox;
  bool m_bShowFaces;
  bool m_bShowEdges;
  bool m_bShowCuts;
  bool m_bShowViolationPolytopes;
  bool m_bShowParticles;
  bool m_bClipping;
  bool m_bShowClippingPlane;
  bool m_bSyncedClipping;
  bool m_bSurfacesOnly;
  bool m_bColorByQuality;

  bool m_bOpenGLError;


  std::vector<bool> m_bMaterialFaceLock;
  std::vector<bool> m_bMaterialCellLock;

  float m_shrinkscale;
  float m_4fvBBoxColor[4];
  float m_4fvCutsColor[4];
  float m_4fvParticleColor[4];
  float m_4fvClippingPlane[4];

  GLuint m_cutVBO;
  GLuint m_violVBO;
  GLuint m_meshVBO[3];

  GLuint m_meshVertexCount;
  GLuint m_cutVertexCount;
  GLuint m_violVertexCount;

  std::vector<Cleaver::vec3> *particles;


  void drawOTCell(Cleaver::OTCell *node);
  void drawTree();
  void drawFaces();
  void drawEdges();
  void drawCuts();
  void drawClippingPlane();
  void drawBox(const Cleaver::BoundingBox &box);

  // draw violation regions around vertices
  void drawViolationPolytopesForVertices();
  void drawViolationPolytopeForVertex(int v);
  //void drawSafetyPolytopes();

  // experimental
  void initializeSSAO();
  void drawFacesWithSSAO();

  //void printModelViewProjection();
  void dumpSVGImage(const std::string &filename);


  // worker functions
  void setup_vbos();
  void update_vbos();
  void build_bkgrnd_vbos();
  void build_output_vbos();

  // star adjacency visualization calls
  void drawVertexStar(int v);
  void drawEdgeStar(int e);
  void drawFaceStar(int f);

  void drawParticles();



  GLenum program;
  GLenum vertex_shader;
  GLenum fragment_shader;

  protected:

  void initializeOptions();
  void initializeCamera();
  void initializeShaders();
  void initializeGL();
  void paintGL();
  void resizeGL(int width, int height);
  void mousePressEvent(QMouseEvent *event);
  void mouseReleaseEvent(QMouseEvent *event);
  void mouseMoveEvent(QMouseEvent *event);
  void keyPressEvent(QKeyEvent *event);
  void keyReleaseEvent(QKeyEvent *event);
  void wheelEvent(QWheelEvent *event);
  void closeEvent(QCloseEvent *event);
};

#endif // MESHWINDOW_H
