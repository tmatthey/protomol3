/*  -*- c++ -*-  */
/*
 *
 * Time your algorithms!
 *
 * 0     5     10    15    20    25    30    35    40    45    50 seconds
 * |-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
 *       ^        ^        ^           ^                      ^
 *       |        |        |           |                      |
 * |<--------------------------->|     |     |<--------->|    |
 * start |        |        |   stop    |     start     stop   |
 *       |        |        |           |                      |
 *       | getTime()=>12s  |    getTime()=>25s         getTime()=>35s
 *       |                 |
 *       |                 |
 *     lap();            lap();
 * getLapTime()=>5s  getLapTime()=>15s
 *
 *
 * Copyright (C) 1995, Christoph Streit <streit@iam.unibe.ch>
 *                     Stephan Amann <amann@iam.unibe.ch>
 *                     University of Berne, Switzerland
 *
 * All rights reserved.
 *
 * This software may be freely copied, modified, and redistributed
 * provided that this copyright notice is preserved on all copies.
 *
 * You may not distribute this software, in whole or in part, as part of
 * any commercial product without the express consent of the authors.
 *
 * There is no warranty or other guarantee of fitness of this software
 * for any purpose.  It is provided solely "as is".
 *
 */

#ifndef _TIMER_H
#define  _TIMER_H
#include <protomol/type/Real.h>
#include <protomol/base/Report.h>

namespace ProtoMol {
  //____________________________________________________________________ TimeRep
  /**
   * Class representation and container of Time
   */
  class TimeRep {
    friend class Timer;
public:
    TimeRep();
    TimeRep(Real realTime, Real userTime, Real sysTime);

    /// Elapsed time
    Real getRealTime() const;
    /// Time spent by this task/thread
    Real getUserTime() const;
    /// Time spent on system operations
    Real getSystemTime() const;
    /// User + System time
    Real getProcessTime() const; // User & Sys

    TimeRep operator+(const TimeRep &time) const;
    TimeRep operator-(const TimeRep &time) const;
    TimeRep &operator+=(const TimeRep &time);
    TimeRep &operator-=(const TimeRep &time);

    friend Report::MyStreamer &operator<<(Report::MyStreamer &os,
                                          const TimeRep &time);

private:
    Real myRealTime;
    Real myUserTime;
    Real mySystemTime;

private:
    void set(Real realTime, Real userTime, Real sysTime);
    void reset();
  };

  //____________________________________________________________ INLINE TimerRep

  inline Real TimeRep::getRealTime() const {
    return myRealTime;
  }

  inline Real TimeRep::getUserTime() const {
    return myUserTime;
  }

  inline Real TimeRep::getSystemTime() const {
    return mySystemTime;
  }

  inline Real TimeRep::getProcessTime() const {
    return myUserTime + mySystemTime;
  }

  //______________________________________________________________________ Timer

  class Timer {
public:
    Timer();

    void start();
    void stop();
    TimeRep lap();
    void reset();

    TimeRep getTime() const;
    TimeRep getLapTime() const;
    TimeRep getActualTime() const;

    static TimeRep getCurrentTime();

    Timer &operator+=(const TimeRep &time);
    Timer &operator-=(const TimeRep &time);

    friend Report::MyStreamer &operator<<(Report::MyStreamer &os,
                                          const Timer &timer);

private:
    bool myRunningFlag;
    TimeRep myStartTime;
    TimeRep myTotalTime;
    TimeRep myLastLapTime;
    TimeRep myLapTime;
  };

  
  Report::MyStreamer &operator<<(Report::MyStreamer &os, const Timer &timer);
}
#endif // _TIMER_H
