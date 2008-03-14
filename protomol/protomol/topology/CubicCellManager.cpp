#include <protomol/topology/CubicCellManager.h>

#include <protomol/base/Exception.h>

using namespace std;
using namespace ProtoMol;

//____ CubicCellManager

const string CubicCellManager::keyword("Cubic");

CubicCellManager::CubicCellManager(Real r) :
  myCellSize(r), myRealCellSize(Vector3D(r, r, r)),
  myRealRCellSize(Vector3D(1.0 / r, 1.0 / r, 1.0 / r)) {}

void CubicCellManager::setCellSize(Real newSize) {
  myCellSize = newSize;
  myRealCellSize = Vector3D(newSize, newSize, newSize);
  myRealRCellSize = Vector3D(1.0 / newSize, 1.0 / newSize, 1.0 / newSize);
}

void CubicCellManager::initialize(CellListStructure &cellList,
                                  const Vector3D &min, const Vector3D &max,
                                  bool pbc) const {
  if (pbc) {
    Vector3D d(max - min);
    myRealCellSize = Vector3D(d.x / std::max(1.0, floor(d.x / myCellSize)),
      d.y / std::max(1.0, floor(d.y / myCellSize)),
      d.z / std::max(1.0, floor(d.z / myCellSize)));
    myRealRCellSize.x = 1.0 / myRealCellSize.x;
    myRealRCellSize.y = 1.0 / myRealCellSize.y;
    myRealRCellSize.z = 1.0 / myRealCellSize.z;
    cellList.initialize(max - min, myRealCellSize);
  } else
    cellList.initialize(max - min + myRealCellSize, myRealCellSize);
}

void CubicCellManager::updateCache(CellListStructure &cellList) const {
  cellList.updateCache();
}

void CubicCellManager::getParameters(vector<Parameter> &parameters) const {
  parameters.push_back
  (Parameter("cellSize", Value(myCellSize, ConstraintValueType::Positive()),
      Text("For Periodic BC this must be < least cell basis vector."
           "  Typically 1/2 of the least cutoff value")));
}

CubicCellManager CubicCellManager::make(vector<Value> values) {
  Real r = 1.0;

  if (!values[0].get(r) || r <= 0.0)
    THROW(keyword + " cellmanager: cutoff > 0 (" + values[0].getString() +
      ").");

  return CubicCellManager(r);
}

