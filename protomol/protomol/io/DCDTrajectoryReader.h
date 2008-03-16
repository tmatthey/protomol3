/*  -*- c++ -*-  */
#ifndef DCDTRAJECTORYREADER_H
#define DCDTRAJECTORYREADER_H

#include <protomol/io/Reader.h>
#include <protomol/type/XYZ.h>
#include <protomol/type/TypeSelection.h>

#include <vector>

namespace ProtoMol {
  //____DCDTrajectoryReader

  /**
   * Reads a DCD trajectory file, frame by frame. Automatic endianess
   * detection.
   */
  class DCDTrajectoryReader : public Reader {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Typedef
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    typedef TypeSelection::Int<4>::type int32;
    typedef TypeSelection::Float<4>::type float4;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors (both default here), assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    DCDTrajectoryReader();
    explicit DCDTrajectoryReader(const std::string &filename);
    virtual ~DCDTrajectoryReader();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Reader
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool tryFormat();
    virtual bool read();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class DCDTrajectory
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

    bool read(std::vector<XYZ> &coords);
    void doRead(std::vector<XYZ> &coords);

    std::vector<XYZ> *orphanXYZ();

  private:
    void fortranRead(char *data, unsigned int size,
                     const std::string &err = "");
    char *fortranReadX(char *data, unsigned int &size,
                      const std::string &err = "");

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Friends
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    friend DCDTrajectoryReader &operator>>(DCDTrajectoryReader &reader,
                                           std::vector<XYZ> &xyz);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    std::vector<XYZ> *xyz;
    bool swap;
  };
}
#endif /* DCDTRAJECTORYREADER_H */

