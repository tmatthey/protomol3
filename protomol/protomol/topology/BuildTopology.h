#ifndef BUILD_TOPOLOGY_H
#define BUILD_TOPOLOGY_H

#include <protomol/topology/ExclusionType.h>

namespace ProtoMol {
  class GenericTopology;
  class PSF;
  class PAR;

  void buildExclusionTable(GenericTopology *topo,
                           const ExclusionType &exclusionType);
  void buildTopology(GenericTopology *topo, const PSF &psf,
                     const PAR &par, bool dihedralMultPSF);
  void buildMoleculeTable(GenericTopology *topo);
}

#endif // BUILD_TOPOLOGY_H
