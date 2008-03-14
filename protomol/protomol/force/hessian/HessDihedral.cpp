#include <protomol/force/hessian/HessDihedral.h>
#include <protomol/type/Vector3DBlock.h>

using namespace ProtoMol;
//constructors
HessDihedral::HessDihedral() {}

HessDihedral::HessDihedral(const Torsion &currTorsion, //torsion data
                           const Vector3D &a1, const Vector3D &a2,  //positions
                           const Vector3D &a3, const Vector3D &a4) {
  evaluate(currTorsion, a1, a2, a3, a4);
  // call the evaluate function
}

void HessDihedral::evaluate(const Torsion &currTorsion, const Vector3D &a1,
                            const Vector3D &a2, const Vector3D &a3,
                            const Vector3D &a4) {
  //actual positions
  double a[9] = {
    a1[0], a2[0], a3[0], a1[1], a2[1], a3[1], a1[2], a2[2], a3[2]
  };
  //determinant of position matrix
  double dta = a[0] * (a[8] * a[4] - a[7] * a[5]) - a[3] *
    (a[8] * a[1] - a[7] * a[2]) + a[6] * (a[5] * a[1] - a[4] * a[2]);
  //matrix will be singular if forst 3 atoms already in a plane containing one
  // of the axes probability is very small so ignore in this case.
  if (dta == 0.0) {
    for (int i = 0; i < 144; i++) hessD[i] = 0.0;

    //clear Hessian
    return;
  }
  //Vector3D rxy // Vector from atom a to atom b
  Vector3D r12(a2 - a1);
  Vector3D r23v(a3 - a2);

  // Cross product of r12 and r23v, represents the plane shared by these two
  // vectors
  double cosPhi = -r12.dot(r23v) / (r12.norm() * r23v.norm());
  double sinPhi = (r12.cross(r23v)).norm() / (r12.norm() * r23v.norm());
  //  sin(acos(cosPhi));

  //double sinPhi= sin(acos(cosPhi));
  double nsa1 = a1.normSquared(); double nsa2 = a2.normSquared();
  double nsa3 = a3.normSquared();
  //
  double x, y, z;
  z = (nsa3 - nsa2 - r23v.norm() * r23v.norm()) / (-2.0 * r23v.norm());
  x =
    (nsa1 - nsa2 - r12.norm() * r12.norm() + 2.0 * z * r12.norm() *
     cosPhi) / (-2.0 * r12.norm() * sinPhi);
  y = sqrt(nsa2 - x * x - z * z);
  //target positions
  double opC[9] = {
    x - r12.norm() * sinPhi, x, x, y, y, y, z - r12.norm() * cosPhi, z, z -
    r23v.norm()
  };
  //find determinant of opC
  double dtopC = opC[0] *
                 (opC[8] * opC[4] - opC[7] *
                  opC[5]) - opC[3] *
                 (opC[8] * opC[1] - opC[7] *
                  opC[2]) + opC[6] * (opC[5] * opC[1] - opC[4] * opC[2]);
  //fix sign of y if affine transformation not pure rotation (reflection),
  // det must be 1.0 not -1.0
  if (dtopC * dta < 0.0) {   //test should be dtopC/dta == -1, but not safe
    opC[3] *= -1.0; opC[4] *= -1.0; opC[5] *= -1.0;
  }
  //
  //|a11 a12 a13|-1        |  a33a22-a32a23  -(a33a12-a32a13)   a23a12-a22a13 |
  //|a21 a22 a23|=  1/DET *| -(a33a21-a31a23)  a33a11-a31a13 -(a23a11-a21a13) |
  //|a31 a32 a33|          |  a32a21-a31a22  -(a32a11-a31a12)   a22a11-a21a12 |
  //inverse of a
  double ia[9] = {
    a[8] * a[4] - a[7] * a[5],
    -(a[8] * a[1] - a[7] * a[2]), a[5] * a[1] - a[4] * a[2],
    -(a[8] * a[3] - a[6] * a[5]), a[8] * a[0] - a[6] * a[2],
    -(a[5] * a[0] - a[3] * a[2]), a[7] * a[3] - a[6] * a[4],
    -(a[7] * a[0] - a[6] * a[1]), a[4] * a[0] - a[3] * a[1]
  };
  for (int i = 0; i < 9; i++) ia[i] /= dta;

  //rotation matrix
  double aRot[9];
  for (int i = 0; i < 9; i++) aRot[i] = 0.0;

  for (int i = 0; i < 9; i++)
    for (int j = 0; j < 3; j++) aRot[i] +=
        opC[(i / 3) * 3 + j] * ia[i % 3 + j * 3];

  //fourth body position
  double in4[3] = {
    a4[0], a4[1], a4[2]
  };
  //rotate it
  double out4[3] = {
    0.0, 0.0, 0.0
  };
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) out4[i] += aRot[i * 3 + j] * in4[j];

  double dgd[12] = {
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
  };
  // Calculate phi
  //#######################
  double g;
  g = atan2((out4[1] - opC[5]), -(out4[0] - opC[2]));

  //multiplicity
  double fact1 = 0.0;
  double fact2 = 0.0;
  for (int dm = 0; dm < currTorsion.multiplicity; dm++)  //Dihedrals
    if (currTorsion.periodicity[dm] > 0) {
      fact1 += -currTorsion.forceConstant[dm] *
               sin(currTorsion.periodicity[dm] * g +
                   currTorsion.phaseShift[dm]) *
               currTorsion.periodicity[dm];
      fact2 += -currTorsion.forceConstant[dm] *
               cos(currTorsion.periodicity[dm] * g +
                   currTorsion.phaseShift[dm]) *
               currTorsion.periodicity[dm] * currTorsion.periodicity[dm];
    } else {    //Impropers
      Real diff = g - currTorsion.phaseShift[dm];
      if (diff < -M_PI)
        diff += 2 * M_PI;
      else if (diff > M_PI)
        diff -= 2 * M_PI;
      fact1 += 2.0 * currTorsion.forceConstant[dm] * diff;
      fact2 += 2.0 * currTorsion.forceConstant[dm];
    }

  //Setup variables
  Real x21, x21_2, x43, x43_2, x43_3, x43_4, x43_5, y43, y43_2, y43_3, z21,
       z21_2, z43, z43_2, r23, r23_2;
  x21 = opC[1] - opC[0]; x43 = out4[0] - opC[2]; y43 = out4[1] - opC[5];
  z21 = opC[7] - opC[6]; z43 = out4[2] - opC[8]; r23 = r23v.norm();
  x21_2 = x21 * x21; x43_2 = x43 * x43; y43_2 = y43 * y43;  r23_2 = r23 * r23;
  z21_2 = z21 * z21; z43_2 = z43 * z43;
  x43_3 = x43_2 * x43; y43_3 = y43_2 * y43; x43_4 = x43_2 * x43_2; x43_5 =
    x43_2 * x43_3;
  //a1 -force_x=0
  //   -force_y=-1/x21
  dgd[1] = -1 / x21;
  //a4 -force_x=y43/(x43_2*(1+y43_2/x43_2))
  //   -force_y=-1/(x43*(1+y43_2/x43_2))
  dgd[9] = y43 / (x43_2 * (1 + y43_2 / x43_2));
  dgd[10] = -1 / (x43 * (1 + y43_2 / x43_2));
  //a2 -force_x=-y43*z43/(x43_2*r23*(1+y43_2/x43_2))
  //   -force_y=(1-z21/r23)/x21+z43/(r23*x43*(1+y43_2/x43_2))
  dgd[3] = -y43 * z43 / (x43_2 * r23 * (1 + y43_2 / x43_2));
  dgd[4] = (1 - z21 / r23) / x21 + z43 / (r23 * x43 * (1 + y43_2 / x43_2));
  //a3 -force_x=y43*(-1+z43/r23)/(x43_2*(1+y43_2/x43_2))
  //   -force_y=z21/(r23*x21)-(-1+z43/r23)/(x43*(1+y43_2/x43_2))
  dgd[6] = y43 * (-1 + z43 / r23) / (x43_2 * (1 + y43_2 / x43_2));
  dgd[7] = z21 / (r23 * x21) - (-1 + z43 / r23) / (x43 * (1 + y43_2 / x43_2));
  //
  //Hessian
  //dfdxdy := -k*cos(n*g(x, y, z)+ps)*n^2*(diff(g(x, y, z), y))*
  //          (diff(g(x, y, z), x))
  //				-k*sin(n*g(x, y, z)+ps)*n*(diff(g(x, y, z), y, x))
  //create output hessian in x_1 y_1 z_1
  //                         y_1
  //                         z_1
  //format
  for (int i = 0; i < 144; i++) hessD[i] = 0.0;

  //
  //common factors
  double sqTerm1 = (1 + y43_2 / x43_2) * (1 + y43_2 / x43_2);
  double sqTerm2 = (-1 + z43 / r23) * (-1 + z43 / r23);
  //#d2x1y1
  hessD[0 * 12 + 1] = hessD[1 * 12 + 0] = -1 / x21_2;
  //d2x1y2
  hessD[0 * 12 + 3 + 1] = hessD[(3 + 1) * 12 + 0] = (1 - z21 / r23) / x21_2;
  //d2x1y3
  hessD[0 * 12 + 6 + 1] = hessD[(6 + 1) * 12 + 0] = z21 / (r23 * x21_2);
  //**d2x2
  hessD[3 * 12 + 3] = -2 * y43 * z43_2 / (x43_3 * r23_2 * (1 + y43_2 / x43_2)) -
    y43 / (x43 * r23_2 * (1 + y43_2 / x43_2)) + 2 * y43_3 * z43_2 /
    (x43_5 * r23_2 * sqTerm1); //changed 2*y43/
  //**d2x3
  hessD[6 * 12 + 6] = -2 * y43 * sqTerm2 / (x43_3 * (1 + y43_2 / x43_2)) -
    y43 / (x43 * r23_2 * (1 + y43_2 / x43_2)) + 2 * y43_3 * sqTerm2 /
    (x43_5 * sqTerm1);               //changed from 2*y43/(
  //#d2x4
  hessD[9 * 12 + 9] = -2 * y43 / (x43_3 * (1 + y43_2 / x43_2)) + 2 * y43_3 /
    (x43_5 * sqTerm1);
  //**d2y2
  hessD[(3 + 1) * 12 + 3 + 1] = y43 / (x43 * r23_2 * (1 + y43_2 / x43_2)) +
    2 * z43_2 * y43 / (r23_2 * x43_3 * sqTerm1); //2*
  //**d2y3
  hessD[(6 + 1) * 12 + 6 + 1] = y43 / (x43 * r23_2 * (1 + y43_2 / x43_2)) +
    2 * sqTerm2 * y43 / (x43_3 * sqTerm1);               //2*
  //d2y4
  hessD[(9 + 1) * 12 + 9 + 1] = 2 * y43 / (x43_3 * sqTerm1);
  //
  //**d2x2x3
  hessD[(3) * 12 + 6] = hessD[(6) * 12 + 3] =
    2 * y43 * z43 * (-1 + z43 / r23) / (x43_3 * r23 * (1 + y43_2 / x43_2)) +
    y43 / (x43 * r23_2 * (1 + y43_2 / x43_2)) - 2 * y43_3 * z43 *
    (-1 + z43 / r23) / (x43_5 * r23 * sqTerm1);  //2*
  //#d2x2x4
  hessD[(3) * 12 + 9] = hessD[(9) * 12 + 3] = 2 * y43 * z43 /
    (x43_3 * r23 * (1 + y43_2 / x43_2)) - 2 * y43_3 * z43 /
    (x43_5 * r23 * sqTerm1);
  //%%%%d2x2y1
  hessD[(3) * 12 + 1] = hessD[(1) * 12 + 3] = (1 - z21 / r23) / x21_2;
  //changed from 1/x21_2
  //%%%%d2x2y2
  hessD[(3) * 12 + 3 + 1] = hessD[(3 + 1) * 12 + 3] =
    -(1 - z21 / r23) * (1 - z21 / r23) / x21_2 + z43_2 /
    (r23_2 * x43_2 * (1 + y43_2 / x43_2)) - 2 * y43_2 * z43_2 /
    (x43_4 * r23_2 * sqTerm1); //changed from -(1-z21/r23)/x21_2+
  //**d2x2y3
  hessD[(3) * 12 + 6 + 1] = hessD[(6 + 1) * 12 + 3] = -z21 *
    (1 - z21 / r23) / (r23 * x21_2) - (-1 + z43 / r23) * z43 /
    (x43_2 * r23 * (1 + y43_2 / x43_2)) + 2 * y43_2 * z43 *
    (-1 + z43 / r23) / (x43_4 * r23 * sqTerm1);
  //#d2x2y4
  hessD[(3) * 12 + 9 + 1] = hessD[(9 + 1) * 12 + 3] = -z43 /
    (x43_2 * (1 + y43_2 / x43_2) * r23) + 2 * y43_2 * z43 /
    (x43_4 * sqTerm1 * r23);                        //=x4y2
  //
  //#d2x3x4
  hessD[(6) * 12 + 9] = hessD[(9) * 12 + 6] = -2 * y43 * (-1 + z43 / r23) /
    (x43_3 * (1 + y43_2 / x43_2)) + 2 * y43_3 * (-1 + z43 / r23) /
    (x43_5 * sqTerm1);
  //d2x3y2!!!!!!!!!!!
  hessD[(6) * 12 + 3 + 1] = hessD[(3 + 1) * 12 + 6] =
    -z43 * (-1 + z43 / r23) / (r23 * x43_2 * (1 + y43_2 / x43_2)) + 2 * y43_2 *
    (-1 + z43 / r23) * z43 / (x43_4 * sqTerm1 * r23) - z21 * (1 - z21 / r23) /
    (x21_2 * r23);
  //d2x3y3!!!!!!!!!!!
  hessD[(6) * 12 + 6 + 1] = hessD[(6 + 1) * 12 + 6] =
    sqTerm2 / (x43_2 * (1 + y43_2 / x43_2)) - 2 * y43_2 * sqTerm2 /
    (x43_4 * sqTerm1) - z21_2 / (r23_2 * x21_2);
  //#d2x3y4
  hessD[(6) * 12 + 9 + 1] = hessD[(9 + 1) * 12 + 6] =
    (-1 + z43 / r23) / (x43_2 * (1 + y43_2 / x43_2)) - 2 * y43_2 *
    (-1 + z43 / r23) / (x43_4 * sqTerm1);           //=x4y3
  //##
  //#d2x4y2
  hessD[(9) * 12 + 3 + 1] = hessD[(3 + 1) * 12 + 9] = -z43 /
    (r23 * x43_2 * (1 + y43_2 / x43_2)) + 2 * y43_2 * z43 /
    (x43_4 * sqTerm1 * r23);                        //=x2y4
  //#d2x4y3
  hessD[(9) * 12 + 6 + 1] = hessD[(6 + 1) * 12 + 9] = hessD[(6) * 12 + 9 + 1];
  //d2x3y4;
  //(-1+z43/r23)/(x43^2*(1+y43^2/x43^2))-2*y43^2*(-1+z43/r23)/
  //  (x43^4*(1+y43^2/x43^2)^2)
  //=x3y4
  //#d2x4y4
  hessD[(9) * 12 + 9 + 1] = hessD[(9 + 1) * 12 + 9] = 1 /
    (x43_2 * (1 + y43_2 / x43_2)) - 2 * y43_2 / (x43_4 * sqTerm1);
  //
  //**d2y2y3
  hessD[(3 + 1) * 12 + 6 + 1] = hessD[(6 + 1) * 12 + 3 + 1] =
    -y43 / (x43 * r23_2 * (1 + y43_2 / x43_2)) - 2 * z43 * y43 *
    (-1 + z43 / r23) / (r23 * x43_3 * sqTerm1);   //2*
  //#d2y2y4
  hessD[(3 + 1) * 12 + 9 + 1] =
    hessD[(9 + 1) * 12 + 3 + 1] = -2 * z43 * y43 / (r23 * x43_3 * sqTerm1);
  //
  //#d2y3y4
  hessD[(6 + 1) * 12 + 9 + 1] = hessD[(9 + 1) * 12 + 6 + 1] = 2 *
    (-1 + z43 / r23) * y43 / (x43_3 * sqTerm1);
  //
  //**d2x2z3
  double temphd = y43 / (x43_2 * r23 * (1 + y43_2 / x43_2));      //same but +/-
  double temphdn = y43 * (1 - z43 / r23) / (x43_2 * r23 * (1 + y43_2 / x43_2));
  hessD[(3) * 12 + 6 + 2] = hessD[(6 + 2) * 12 + 3] = temphdn;   //temphd;
  //#d2x2z4
  hessD[(3) * 12 + 9 + 2] = hessD[(9 + 2) * 12 + 3] = -temphd;   //-d2x2z3;
  //-y43/(x43^2*r23*(1+y43^2/x43^2))
  //**d2x3z3
  hessD[(6) * 12 + 6 + 2] = hessD[(6 + 2) * 12 + 6] = -temphdn;   //-temphd;
  //-d2x2z3; //-y43/(x43^2*r23*(1+y43^2/x43^2))
  //#d2x3z4
  hessD[(6) * 12 + 9 + 2] = hessD[(9 + 2) * 12 + 6] = temphd;   //d2x2z3;
  //y43/(x43^2*r23*(1+y43^2/x43^2))
  //
  //d2y2z1
  temphd = 1 / (r23 * x21);
  hessD[(3 + 1) * 12 + 2] = hessD[(2) * 12 + 3 + 1] = temphd;   //1/(r23*x21);
  //**Temp for y2z2 and -y3z2
  double temphd3 =
    -(1 / r23 - z21 / r23_2) / x21 - z43 / (r23_2 * x43 * (1 + y43_2 / x43_2));
  //** Tempf for y3z3 and -y2z3
  double temphd4 = z21 / (r23_2 * x21) -
    (-1 / r23 + z43 / r23_2) / (x43 * (1 + y43_2 / x43_2));
  //**d2y2z2
  hessD[(3 + 1) * 12 + 3 + 2] = hessD[(3 + 2) * 12 + 3 + 1] = temphd3;
  //-temphd; //-d2y2z1; //-1/(r23*x21)
  //#d2y2z4
  double temphd2 = 1 / (r23 * x43 * (1 + y43_2 / x43_2));
  hessD[(3 + 1) * 12 + 9 + 2] = hessD[(9 + 2) * 12 + 3 + 1] = temphd2;
  //**d2y2z3
  hessD[(3 + 1) * 12 + 6 + 2] = hessD[(6 + 2) * 12 + 3 + 1] = -temphd4;
  //-temphd2; //-d2y2z4; //-1/(r23*x43*(1+y43^2/x43^2))
  //
  //d2y3z1
  hessD[(6 + 1) * 12 + 2] = hessD[(2) * 12 + 6 + 1] = -temphd;   //-d2y2z1;
  //-1/(r23*x21)
  //**d2y3z2
  hessD[(6 + 1) * 12 + 3 + 2] = hessD[(3 + 2) * 12 + 6 + 1] = -temphd3;
  //temphd; //d2y2z1; // 1/(r23*x21)
  //**d2y3z3
  hessD[(6 + 1) * 12 + 6 + 2] = hessD[(6 + 2) * 12 + 6 + 1] = temphd4;
  //temphd2; //d2y2z4; //1/(r23*x43*(1+y43^2/x43^2))
  //#d2y3z4
  hessD[(6 + 1) * 12 + 9 + 2] = hessD[(9 + 2) * 12 + 6 + 1] = -temphd2;
  //-d2y2z4; //-1/(r23*x43*(1+y43^2/x43^2))
  //d2x3y1!!!!!!!!!!
  hessD[(6) * 12 + 1] = hessD[(1) * 12 + 6] = z21 / (x21_2 * r23);
  ///////
  int eye, eyeh, jay, kay;
  for (int i = 0; i < 4; i++)
    //####for test! should be j=i; so 1/2 matrix
    for (int j = 0; j < 4; j++) {
      //create the i,j hessian matrix
      eye = i * 3; eyeh = eye * 12; jay = j * 3;
      for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++) {
          kay = eyeh + jay + k * 12 + l;        // (k/3)*12 + jay + k%3;
          hessD[kay] = hessD[kay] * fact1 +
            dgd[eye + k] * dgd[jay + l] * fact2;
        }

      for (int k = 0; k < 3; k++) {
        double v[3] = {
          hessD[eyeh + jay + k], hessD[eyeh + jay + k + 12],
          hessD[eyeh + jay + k + 24]
        };
        rotateV3D(aRot, v);
        hessD[eyeh + jay + k] = v[0]; hessD[eyeh + jay + k + 12] = v[1];
        hessD[eyeh + jay + k + 24] = v[2];
      }

      for (int k = 0; k < 3; k++)
        rotateV3D(aRot, &hessD[eyeh + jay + 12 * k]);
    }

}

//Use aRot to rotate the vector back into real space
double *HessDihedral::rotateV3D(double *aRot, double *mf) {
  double out[3] = {
    0.0, 0.0, 0.0
  };
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) out[i] += aRot[i + j * 3] * mf[j];

  for (int i = 0; i < 3; i++) mf[i] = out[i];

  return mf;
}
