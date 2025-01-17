/*******************************************************************\

              Copyright (C) 2003 Joseph Coffland

    This program is free software; you can redistribute it and/or
     modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
        of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
             GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
     along with this program; if not, write to the Free Software
      Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
                           02111-1307, USA.

            For information regarding this software email
                   jcofflan@users.sourceforge.net

\*******************************************************************/
#ifndef STRING_H
#define STRING_H

#include <string>
#include <protomol/type/Real.h>
namespace ProtoMol {

  /**
   * Used for convenient conversion of basic data types to std::string.
   */
  class String : public std::string {
  public:
    // Pass through constructors
    /// See std::string
    String(const char *s, size_type n) : std::string(s, n) {}
  
    /// See std::string
    String(const char *s) : std::string(s) {}
  
    /// See std::string
    String(std::string &s) : std::string(s) {}
  
    /// See std::string
    String(std::string &s, size_type pos, size_type n = npos) :
      std::string(s, pos, n) {}
  
    /// See std::string
    String(size_type n, char c) : std::string(n, c) {}
  
    // Conversion constructors
    /** 
     * Convert an int to a character string.
     * 
     * @param x The int to convert.
     */
    String(const int x);

    /** 
     * Convert an unsigned int to a string.
     * 
     * @param x The unsigned int to convert.
     */
    String(const unsigned x);

    /** 
     * Convert a long to a string.
     * 
     * @param x The long to convert.
     */
    String(const long x);

    /** 
     * Convert an unsigned long to a string.
     * 
     * @param x The unsigned long to convert.
     */
    String(const unsigned long x);

    /** 
     * Convert a double to a string.  
     * 
     * @param x The double to convert.
     */
    String(const Real x);

    static unsigned char parseUByte(const std::string s);
    static char parseByte(const std::string s);
    static unsigned short parseUShort(const std::string s);
    static short parseShort(const std::string s);
    static unsigned int parseUInteger(const std::string s);
    static int parseInteger(const std::string s);
    static double parseDouble(const std::string s);
    static bool parseBool(const std::string s);  
    static std::string trim(const std::string s);
    static std::string toUpper(const std::string s);
    static std::string toLower(const std::string s);
  };
}
#endif // STRING_H
