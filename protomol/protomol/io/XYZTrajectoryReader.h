/*  -*- c++ -*-  */
#ifndef XYZTRAJECTORYREADER_H
#define XYZTRAJECTORYREADER_H

#include <protomol/io/XYZReader.h>
#include <protomol/type/XYZ.h>

namespace ProtoMol {
  //____XYZTrajectoryReader

  /**
   * Reads a XYY trajectory files (ASCII).
   */
  class XYZTrajectoryReader : public XYZReader {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors (both default here), assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    XYZTrajectoryReader();
    explicit XYZTrajectoryReader(const std::string &filename);
    virtual ~XYZTrajectoryReader();

    virtual bool tryFormat();
    virtual bool read();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class XYZTrajectoryReader
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    bool read(std::vector<XYZ> &xyz);
    void doRead(std::vector<XYZ> &xyz);

    std::vector<XYZ> *orphanXYZ();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Friends
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    friend XYZTrajectoryReader &operator>>(XYZTrajectoryReader &reader,
                                           std::vector<XYZ> &xyz);

  private:
    std::vector<XYZ> *xyz;
  };
}
#endif /* XYZTRAJECTORYREADER_H */
